#include "microphysics.hpp"
#include "../parameter_input.hpp"
#include "../utils/utils.hpp"
#include "../mesh/mesh.hpp"

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

  std::vector<std::string> str;

  // read eps, size should be NVAPOR
  if (NVAPOR > 0)
    SplitString(pin->GetOrAddString("microphysics", "eps", ""), str);
  for (int n = 0; n < ITR; ++n) eps_[n] = 1.;
  for (int n = 0; n < std::min(NVAPOR, (int)str.size()); ++n) {
    eps_[1+n] = atof(str[n].c_str());
    eps_[1+NVAPOR+n] = atof(str[n].c_str());
  }

  // read rcv, size should be NVAPOR
  if (NVAPOR > 0)
    SplitString(pin->GetOrAddString("microphysics", "rcv", ""), str);
  for (int n = 0; n < ITR; ++n) rcv_[n] = 1.;
  for (int n = 0; n < std::min(NVAPOR, (int)str.size()); ++n) {
    rcv_[1+n] = atof(str[n].c_str());
    rcv_[1+NVAPOR+n] = atof(str[n].c_str());
  }

  // read latent, size should be NVAPOR, unit is [J/kg]
  if (NVAPOR > 0)
    SplitString(pin->GetOrAddString("microphysics", "latent", ""), str);
  for (int n = 0; n < ITR; ++n) latent_[n] = 0.;
  for (int n = 0; n < std::min(NVAPOR, (int)str.size()); ++n) {
    latent_[1+n] = atof(str[n].c_str());
    latent_[1+NVAPOR+n] = atof(str[n].c_str());
  }

  // terminal velocity
  tv_ = pin->GetOrAddReal("microphysics", "terminalv", 0.);

  // auto conversion time
  autoc_ = pin->GetOrAddReal("microphysics", "autoc", 0.);

  // tiny number
  tiny_number_ = pin->GetOrAddReal("microphysics", "tiny_number", 1.E-20);
}

Microphysics::~Microphysics()
{}

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
