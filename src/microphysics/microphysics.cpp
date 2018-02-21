#include <algorithm>  // fill

#include "microphysics.hpp"
#include "../parameter_input.hpp"
#include "../utils/utils.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../math_funcs.hpp"

Microphysics::Microphysics(MeshBlock *pmb, ParameterInput *pin)
{
  pmy_block_ = pmb;
  if (NVAPOR > 0) {
    Rd_ = pin->GetReal("microphysics", "Rd");
    ncycle = pin->GetOrAddInteger("microphysics", "ncycle", 1);
  } else {
    Rd_ = 0.;
    ncycle = 0.;
  }

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
    latent_[n+NVAPOR] = - beta_[n]*Rd_/eps_[n]*t3_[n];
    rcp_[n+NVAPOR] = rcp_[n] + beta_[n+NVAPOR]/eps_[n]*(1. - 1./gamma);
    rcv_[n+NVAPOR] = rcp_[n+NVAPOR]*gamma;
  }

  // auto-conversion time
  autoc_ = pin->GetOrAddReal("microphysics", "autoc", 2000.);

  // evaporation rate
  evapr_ = pin->GetOrAddReal("microphysics", "evapr", 10.);

  // terminal velocity
  termv_ = pin->GetOrAddReal("microphysics", "termv", -10.);

  // allocate temperature
  int ncells1 = pmb->block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  if (pmb->block_size.nx2 > 1) ncells2 = pmb->block_size.nx2 + 2*(NGHOST);
  if (pmb->block_size.nx3 > 1) ncells3 = pmb->block_size.nx3 + 2*(NGHOST);

  T.NewAthenaArray(ncells3, ncells2, ncells1);
  P.NewAthenaArray(ncells3, ncells2, ncells1);
  recondense_.NewAthenaArray(NVAPOR,ncells3,ncells2,ncells1);
  std::fill(recondense_.data(), recondense_.data() + recondense_.GetSize(), true);

  /* debug
  std::cout << "eps_ = ";
  for (int n = 0; n < ITR; ++n)
    std::cout << eps_[n] << " ";
  std::cout << std::endl;

  std::cout << "beta_ = ";
  for (int n = 0; n < ITR; ++n)
    std::cout << beta_[n] << " ";
  std::cout << std::endl;

  std::cout << "rcp_ = ";
  for (int n = 0; n < ITR; ++n)
    std::cout << rcp_[n] << " ";
  std::cout << std::endl;

  std::cout << "rcv_ = ";
  for (int n = 0; n < ITR; ++n)
    std::cout << rcv_[n] << " ";
  std::cout << std::endl;

  std::cout << "latent_ = ";
  for (int n = 0; n < ITR; ++n)
    std::cout << latent_[n] << " ";
  std::cout << std::endl;*/
}

Microphysics::~Microphysics()
{
  T.DeleteAthenaArray();
  P.DeleteAthenaArray();
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
          LE += latent_[n]*u(n,k,j,i);
        T(k,j,i) = (u(IEN,k,j,i) - KE - LE)/cv;

        P(k,j,i) = 0.;
        for (int n = 0; n < 1 + NVAPOR; ++n)
          P(k,j,i) += Rd_*T(k,j,i)*u(n,k,j,i)/eps_[n];
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

Real Microphysics::Tempv(AthenaArray<Real> const& w, int i, int j, int k) const
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
  return Tempv(w,i,j,k)/feps;
}

Real Microphysics::Theta(Real p0, AthenaArray<Real> const& w, int i, int j, int k) const
{
  Real chi = Chi(w,i,j,k);
  Real temp = Temp(w,i,j,k);
  return temp*pow(p0/w(IPR,k,j,i), chi);
}

Real Microphysics::Thetav(Real p0, AthenaArray<Real> const& w, int i, int j, int k) const
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
    LE += latent_[n]*w(n,k,j,i);
  return Cp(w,i,j,k)*Temp(w,i,j,k) + LE + grav*pcoord->x1v(i);
}

// Dropping Precipitation
#if PRECIPITATION_ENABLED
void Hydro::TracerAdvection(int k, int j, int il, int iu, int ivx,
  AthenaArray<Real> const& wl, AthenaArray<Real> const& wr, AthenaArray<Real> &flx)
{
  Real tv = pmy_block->pmicro->GetTerminalVelocity();
  for (int i = il; i <= iu; ++i) {
    if (ivx == IVX) {
      for (int n = ITR; n < ITR + NVAPOR; ++n) {
        flx(n,i) = 0.5*((wl(ivx,i) + tv) + fabs(wl(ivx,i) + tv))*wl(n,i)*wl(IDN,i)
                 + 0.5*((wr(ivx,i) + tv) - fabs(wr(ivx,i) + tv))*wr(n,i)*wr(IDN,i);
      }
      for (int n = ITR + NVAPOR; n < ITR + NTRACER; ++n) {
        flx(n,i) = 0.5*(wl(ivx,i) + fabs(wl(ivx,i)))*wl(n,i)*wl(IDN,i)
                 + 0.5*(wr(ivx,i) - fabs(wr(ivx,i)))*wr(n,i)*wr(IDN,i);
      }
    } else {
      for (int n = ITR; n < ITR + NTRACER; ++n) {
        flx(n,i) = 0.5*(wl(ivx,i) + fabs(wl(ivx,i)))*wl(n,i)*wl(IDN,i)
                 + 0.5*(wr(ivx,i) - fabs(wr(ivx,i)))*wr(n,i)*wr(IDN,i);
      }
    }
  }
}
#endif
