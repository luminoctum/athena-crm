#include <cstdlib>

#include "../athena_arrays.hpp"
#include "../hydro/hydro.hpp"
#include "../math_funcs.hpp"
#include "microphysics.hpp"

void Microphysics::Precipitation(AthenaArray<Real> &u, Real dt)
{
  std::stringstream msg;
  if (dt > autoc_) {
    msg << dt << " " << autoc_ << " ";
    msg << "Microphysics time step smaller than autoconversion time" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  MeshBlock *pmb = pmy_block_;
  Real gamma = pmb->peos->GetGamma();

  //Real max = 0.;
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i)
        for (int n = 0; n < NVAPOR; ++n) {
          Real drho;
          int nc = ICD + n;
          int np = ITR + n;
          if (recondense_(n,k,j,i))
            drho = u(nc,k,j,i);
          else
            drho = u(nc,k,j,i)*dt/autoc_;
          u(nc,k,j,i) -= drho;
          u(np,k,j,i) += drho;
          u(IEN,k,j,i) -= Rd_/(gamma - 1)*drho*rcv_[nc]*temp_(k,j,i) - latent_[nc]*drho;
          //if (u(nc,k,j,i) > max) max = u(nc,k,j,i);
        }
  //std::cout << max << std::endl;
}

// Dropping Precipitation
#if PRECIPITATION_ENABLED
void Hydro::TracerAdvection(int k, int j, int il, int iu, int ivx,
  AthenaArray<Real> const& wl, AthenaArray<Real> const& wr, AthenaArray<Real> &flx)
{
  Real grav = -psrc->GetG1();
  Real dp = pmy_block->pmicro->GetDp();
  //Real tv = -10.;
  for (int i = il; i <= iu; ++i) {
    if (ivx == IVX) {
      Real tv = pmy_block->pmicro->TerminalVelocity(dp, wl(IDN,i) + wr(IDN,i), grav);
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
