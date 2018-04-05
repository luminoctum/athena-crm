#include <algorithm>  // fill

#include "microphysics.hpp"
#include "thermodynamics.hpp"
#include "../parameter_input.hpp"
#include "../utils/utils.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../math_funcs.hpp"

Microphysics::Microphysics(MeshBlock *pmb, ParameterInput *pin)
{
  pmy_block_ = pmb;
  Rd_ = pin->GetOrAddReal("microphysics", "Rd", 0.);
  if (NVAPOR > 0)
    ncycle = pin->GetOrAddInteger("microphysics", "ncycle", 1);
  else
    ncycle = 0;

  Real gamma = pin->GetReal("hydro", "gamma");

  std::vector<std::string> str;

  // read eps, size should be NVAPOR
  if (NVAPOR > 0)
    SplitString(pin->GetOrAddString("microphysics", "eps", ""), str);
  for (int n = 0; n < ITR; ++n) eps_[n] = 1.;
  for (int n = 0; n < std::min(NVAPOR, (int)str.size()); ++n) {
    eps_[1+n] = atof(str[n].c_str());
    eps_[1+NVAPOR+n] = atof(str[n].c_str());
  }

  // read rcp, size should be NVAPOR
  if (NVAPOR > 0)
    SplitString(pin->GetOrAddString("microphysics", "rcp", ""), str);
  rcp_[0] = 1.;
  for (int n = 0; n < NVAPOR; ++n)
    rcp_[1+n] = 1./eps_[1+n];
  for (int n = 0; n < std::min(NVAPOR, (int)str.size()); ++n)
    rcp_[1+n] = atof(str[n].c_str())/eps_[1+n];

  // read beta, size should be 2*NVAPOR
  if (NVAPOR > 0)
    SplitString(pin->GetOrAddString("microphysics", "beta", ""), str);
  for (int n = 0; n < ITR; ++n) beta_[n] = 0.;
  for (int n = 0; n < std::min(2*NVAPOR, (int)str.size()); ++n)
    beta_[1+n] = atof(str[n].c_str());

  // read t3, size should be NVAPOR, unit is K
  if (NVAPOR > 0)
    SplitString(pin->GetOrAddString("microphysics", "t3", ""), str);
  for (int n = 0; n < 1+NVAPOR; ++n) t3_[n] = 0.;
  for (int n = 0; n < std::min(NVAPOR, (int)str.size()); ++n)
    t3_[1+n] = atof(str[n].c_str());

  // read p3, size should be NVAPOR, unit is pa
  if (NVAPOR > 0)
    SplitString(pin->GetOrAddString("microphysics", "p3", ""), str);
  for (int n = 0; n < 1+NVAPOR; ++n) p3_[n] = 0.;
  for (int n = 0; n < std::min(NVAPOR, (int)str.size()); ++n)
    p3_[1+n] = atof(str[n].c_str());

  // calculate rcv
  rcv_[0] = 1.;
  for (int n = 1; n < ICD; ++n)
    rcv_[n] = gamma*rcp_[n] + (1. - gamma)/eps_[n];

  // calculate cc and latent using beta and rcv
  latent_[0] = 0.;
  for (int n = 1; n < ICD; ++n) {
    latent_[n] = 0.;
    latent_[n+NVAPOR] = beta_[n]*Rd_/eps_[n]*t3_[n];
    rcp_[n+NVAPOR] = rcp_[n] + beta_[n+NVAPOR]/eps_[n]*(1. - 1./gamma);
    rcv_[n+NVAPOR] = rcp_[n+NVAPOR]*gamma;
  }

  // auto-conversion time
  autoc_ = pin->GetOrAddReal("microphysics", "autoc", 2000.);

  // mass diffusion coefficient, water(g) - air(g), 25 degree C
  Kw_ = pin->GetOrAddReal("microphysics", "Kw", 2.8E-5);

  // heat diffusion coefficient, air, 300 K
  Kt_ = pin->GetOrAddReal("microphysics", "Kt", 1.9E-5);

  // kinematic viscosity
  Kv_ = pin->GetOrAddReal("microphysics", "Kv", 1.6E-5);

  // particle diameter
  Dp_ = pin->GetOrAddReal("microphysics", "Dp", 1.0E-3);

  // particle density
  rhol_ = pin->GetOrAddReal("microphysics", "rhol", 1.0E3);

  // allocate temperature
  int ncells1 = pmb->block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  if (pmb->block_size.nx2 > 1) ncells2 = pmb->block_size.nx2 + 2*(NGHOST);
  if (pmb->block_size.nx3 > 1) ncells3 = pmb->block_size.nx3 + 2*(NGHOST);

  temp_.NewAthenaArray(ncells3, ncells2, ncells1);
  pres_.NewAthenaArray(ncells3, ncells2, ncells1);
  recondense_.NewAthenaArray(NVAPOR,ncells3,ncells2,ncells1);
  std::fill(recondense_.data(), recondense_.data() + recondense_.GetSize(), true);
}

Microphysics::~Microphysics()
{
  temp_.DeleteAthenaArray();
  pres_.DeleteAthenaArray();
  recondense_.DeleteAthenaArray();
}

void Microphysics::CalculateTP(AthenaArray<Real> const& u)
{
  MeshBlock *pmb = pmy_block_;
  Real gamma = pmy_block_->peos->GetGamma();
  Real cvd = Rd_/(gamma - 1);

  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real KE = 0, LE = 0, rho = 0, cv = 0;
        for (int n = 0; n < ITR; ++n) {
          rho += u(n,k,j,i);
          cv += u(n,k,j,i)*cvd*rcv_[n];
        }
        KE = 0.5*(_sqr(u(IM1,k,j,i)) + _sqr(u(IM2,k,j,i)) + _sqr(u(IM3,k,j,i)))/rho;
        for (int n = ICD; n < ICD + NVAPOR; ++n)
          LE += -latent_[n]*u(n,k,j,i);
        temp_(k,j,i) = (u(IEN,k,j,i) - KE - LE)/cv;

        pres_(k,j,i) = 0.;
        for (int n = 0; n < 1 + NVAPOR; ++n)
          pres_(k,j,i) += Rd_*temp_(k,j,i)*u(n,k,j,i)/eps_[n];
      }
}

Real Microphysics::Chi(AthenaArray<Real> const& w, int i, int j, int k) const
{
  Real gamma = pmy_block_->peos->GetGamma();
  Real feps = 1., fsig = 1.;
  for (int n = ICD; n < ICD + NVAPOR; ++n) {
    feps -= w(n,k,j,i);
    fsig += w(n,k,j,i)*(rcp_[n] - 1.);
  }
  for (int n = 1; n < 1 + NVAPOR; ++n) {
    feps += w(n,k,j,i)*(1./eps_[n] - 1.);
    fsig += w(n,k,j,i)*(rcp_[n] - 1.);
  }
  return (gamma - 1.)/gamma*feps/fsig;
}

Real Microphysics::Cp(AthenaArray<Real> const& w, int i, int j, int k) const
{
  Real gamma = pmy_block_->peos->GetGamma();
  Real fsig = 1.;
  for (int n = 1; n < ITR; ++n)
    fsig += w(n,k,j,i)*(rcp_[n] - 1.);
  return gamma/(gamma - 1.)*Rd_*fsig;
}

Real Microphysics::Cv(AthenaArray<Real> const& w, int i, int j, int k) const
{
  Real gamma = pmy_block_->peos->GetGamma();
  Real fsig = 1.;
  for (int n = 1; n < ITR; ++n)
    fsig += w(n,k,j,i)*(rcv_[n] - 1.);
  return 1./(gamma - 1.)*Rd_*fsig;
}

Real Microphysics::TempV(AthenaArray<Real> const& w, int i, int j, int k) const
{
  return w(IPR,k,j,i)/(w(IDN,k,j,i)*Rd_);
}

Real Microphysics::Temp(AthenaArray<Real> const& w, int i, int j, int k) const
{
  Real feps = 1.;
  for (int n = ICD; n < ICD + NVAPOR; ++n)
    feps -= w(n,k,j,i);
  for (int n = 1; n < 1 + NVAPOR; ++n)
    feps += w(n,k,j,i)*(1./eps_[n] - 1.);
  return TempV(w,i,j,k)/feps;
}

Real Microphysics::ThetaE(Real p0, AthenaArray<Real> const& w, int i, int j, int k) const
{
  Real gamma = pmy_block_->peos->GetGamma();
  Real cpd = Rd_*gamma/(gamma - 1.);
  Real temp = Temp(w,i,j,k);
  Real pres = w(IPR,k,j,i);

  Real qd = 1.;
  for (int n = 1; n < ITR; ++n)
    qd -= w(n,k,j,i);

  Real lv = 0.;
  for (int n = ICD; n < ICD + NVAPOR; ++n)
    lv += GetLatent(n, temp)*w(n-ICD+1,k,j,i);

  Real st = 1.;
  for (int n = ICD; n < ICD + NVAPOR; ++n)
    st += (w(n,k,j,i) + w(n-ICD+1,k,j,i))*(GetCpRatio(n) - 1.);
  Real lv_ov_cpt = lv/(cpd*st*temp);

  Real chi = Rd_/cpd*qd/st;

  Real xv = 0.;
  for (int n = 1; n < ICD; ++n)
    xv += w(n,k,j,i)/qd/eps_[n];

  Real pd = pres/(1. + xv);

  Real rh = 1.;
  for (int n = 1; n < ICD; ++n) {
    Real eta = w(n,k,j,i)/qd/eps_[n];
    Real pv = pres*eta/(1. + xv);
    Real esat = SatVaporPresIdeal(temp/t3_[n], p3_[n], beta_[n], beta_[n+NVAPOR]);
    rh *= pow(pv/esat, -eta*Rd_/(cpd*st));
  }

  return temp*pow(p0/pd, chi)*exp(lv_ov_cpt)*rh;
}

Real Microphysics::Theta(Real p0, AthenaArray<Real> const& w, int i, int j, int k) const
{
  Real chi = Chi(w,i,j,k);
  Real temp = Temp(w,i,j,k);
  return temp*pow(p0/w(IPR,k,j,i), chi);
}

Real Microphysics::ThetaV(Real p0, AthenaArray<Real> const& w, int i, int j, int k) const
{
  Real feps = 1.;
  for (int n = ICD; n < ICD + NVAPOR; ++n)
    feps -= w(n,k,j,i);
  for (int n = 1; n < 1 + NVAPOR; ++n)
    feps += w(n,k,j,i)*(1./eps_[n] - 1.);

  return Theta(p0,w,i,j,k)*feps;
}

Real Microphysics::MSE(Real grav, AthenaArray<Real> const& w, int i, int j, int k) const
{
  Coordinates *pcoord = pmy_block_->pcoord;
  Real LE = 0.;
  for (int n = ICD; n < ICD + NVAPOR; ++n)
    LE += -latent_[n]*w(n,k,j,i);
  return Cp(w,i,j,k)*Temp(w,i,j,k) + LE + grav*pcoord->x1v(i);
}

void Microphysics::Prim2Hydro(Real const prim[], AthenaArray<Real>& w,
  int i, int j, int k) const
{
  // set mass mixing ratio
  Real sum = 1.;
  for (int n = 1; n < ITR; ++n) {
    w(n,k,j,i) = prim[n]*eps_[n];
    sum += prim[n]*(eps_[n] - 1.); 
  }
  for (int n = 1; n < ITR; ++n)
    w(n,k,j,i) /= sum;

  // virtual temperature
  Real feps = 1.;
  for (int n = ICD; n < ICD + NVAPOR; ++n)
    feps -= w(n,k,j,i);
  for (int n = 1; n < 1 + NVAPOR; ++n)
    feps += w(n,k,j,i)*(1./eps_[n] - 1.);
  Real Tv = prim[IDN]*feps;

  w(IPR,k,j,i) = prim[IPR];
  w(IDN,k,j,i) = prim[IPR]/(Rd_*Tv);
}

void Microphysics::Hydro2Prim(Real prim[], AthenaArray<Real> const& w,
  int i, int j, int k) const
{
  Real sum = 1.;
  for (int n = 1; n < ITR; ++n) {
    prim[n] = w(n,k,j,i)/eps_[n];
    sum += w(n,k,j,i)*(1./eps_[n] - 1.);
  }
  for (int n = 1; n < ITR; ++n)
    prim[n] /= sum;

  prim[IPR] = w(IPR,k,j,i);
  prim[IDN] = Temp(w,i,j,k);
}
