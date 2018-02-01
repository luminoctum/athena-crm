//! \file lmars_adiabatic.cpp
//  \brief LMARS based on Chen's 2013 paper

// Athena++ headers
#include "../../hydro.hpp"
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../../eos/eos.hpp"
//#include "../../../math_funcs.hpp"

void Hydro::RiemannSolver(int const k, int const j, int const il, int const iu,
    int const ivx, AthenaArray<Real> const& bx, AthenaArray<Real> &wl,
    AthenaArray<Real> &wr, AthenaArray<Real> &flx)
{
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;

  EquationOfState *peos = pmy_block->peos;
  Real rhobar, pbar, cbar, ubar, hl, hr;
  //Real gamma = pmy_block->peos->GetGamma();
  //Real kappa = 1. - 1./gamma;

  Real rho[1+NVAPOR];
  for (int i = il; i <= iu; ++i) {
    // correction for gamma
    Real gmix = 1.;
    for (int n = 1; n < ITR; ++n)
      gmix += 0.5*(wl(n,i) + wr(n,i))*(peos->GetCvRatio(n) - 1.);
    Real gamma = (peos->GetGamma() - 1. + gmix)/gmix;
    Real kappa = 1. - 1./gamma;

    rhobar = 0.5*(wl(IDN,i) + wr(IDN,i));
    hl = wl(IPR,i)/wl(IDN,i)/kappa + 0.5*(wl(ivx,i)*wl(ivx,i) + wl(ivy,i)*wl(ivy,i) + wl(ivz,i)*wl(ivz,i));
    hr = wr(IPR,i)/wr(IDN,i)/kappa + 0.5*(wr(ivx,i)*wr(ivx,i) + wr(ivy,i)*wr(ivy,i) + wr(ivz,i)*wr(ivz,i));
    cbar = sqrt(0.5*gamma*(wl(IPR,i)+ wr(IPR,i))/rhobar);
    pbar = 0.5*(wl(IPR,i) + wr(IPR,i)) + 0.5*(rhobar*cbar)*(wl(ivx,i) - wr(ivx,i));
    ubar = 0.5*(wl(ivx,i) + wr(ivx,i)) + 0.5/(rhobar*cbar)*(wl(IPR,i) - wr(IPR,i));

    if (ubar > 0.) {
      // volume mixing ratio to mass mixing ratio
      Real qd = 1., sum = 0.;
      for (int n = 1; n < ITR; ++n) {
        rho[n] = wl(n,i)*peos->GetMassRatio(n);
        qd -= wl(n,i);
        sum += rho[n];
      }
      rho[0] = qd;
      sum += qd;

      for (int n = 0; n < ITR; ++n)
        flx(n,i) = ubar*wl(IDN,i)*rho[n]/sum;
      flx(ivx,i) = ubar*wl(IDN,i)*wl(ivx,i) + pbar;
      flx(ivy,i) = ubar*wl(IDN,i)*wl(ivy,i);
      flx(ivz,i) = ubar*wl(IDN,i)*wl(ivz,i);
      flx(IEN,i) = ubar*wl(IDN,i)*hl;
    } else {
      // volume mixing ratio to mass mixing ratio
      Real qd = 1., sum = 0.;
      for (int n = 1; n < ITR; ++n) {
        rho[n] = wr(n,i)*peos->GetMassRatio(n);
        qd -= wr(n,i);
        sum += rho[n];
      }
      rho[0] = qd;
      sum += qd;

      for (int n = 0; n < ITR; ++n)
        flx(n,i) = ubar*wr(IDN,i)*rho[n]/sum;
      flx(ivx,i) = ubar*wr(IDN,i)*wr(ivx,i) + pbar;
      flx(ivy,i) = ubar*wr(IDN,i)*wr(ivy,i);
      flx(ivz,i) = ubar*wr(IDN,i)*wr(ivz,i);
      flx(IEN,i) = ubar*wr(IDN,i)*hr;
    }
  }
}
