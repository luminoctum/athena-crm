//! \file lmars_adiabatic.cpp
//  \brief LMARS based on Chen's 2013 paper

// Athena++ headers
#include "../../hydro.hpp"
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../../microphysics/microphysics.hpp"
#include "../../../eos/eos.hpp"
#include "../../../math_funcs.hpp"

void Hydro::RiemannSolver(int const k, int const j, int const il, int const iu,
    int const ivx, AthenaArray<Real> const& bx, AthenaArray<Real> &wl,
    AthenaArray<Real> &wr, AthenaArray<Real> &flx)
{
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;

  Microphysics *pmicro = pmy_block->pmicro;
  Real rhobar, pbar, cbar, ubar, hl, hr;
  Real gamma = pmy_block->peos->GetGamma();
  // kappa = gamma/(gamma - 1)

  for (int i = il; i <= iu; ++i) {
    // correction for gamma
    // left
    Real fsig = 1., feps = 1., lel = 0., ler = 0.;
    for (int n = ICD; n < ICD + NVAPOR; ++n) {
      lel += -pmicro->GetLatent(n)*wl(n,i);
      fsig += wl(n,i)*(pmicro->GetCvRatio(n) - 1.);
      feps -= wl(n,i);
    }
    for (int n = 1; n < ICD; ++n) {
      fsig += wl(n,i)*(pmicro->GetCvRatio(n) - 1.);
      feps += wl(n,i)*(1./pmicro->GetMassRatio(n) - 1.);
    }
    Real kappal = 1./(gamma - 1.)*fsig/feps;

    // right
    fsig = 1., feps = 1.;
    for (int n = ICD; n < ICD + NVAPOR; ++n) {
      ler += -pmicro->GetLatent(n)*wr(n,i);
      fsig += wr(n,i)*(pmicro->GetCvRatio(n) - 1.);
      feps -= wr(n,i);
    }
    for (int n = 1; n < ICD; ++n) {
      fsig += wr(n,i)*(pmicro->GetCvRatio(n) - 1.);
      feps += wr(n,i)*(1./pmicro->GetMassRatio(n) - 1.);
    }
    Real kappar = 1./(gamma - 1.)*fsig/feps;

    rhobar = 0.5*(wl(IDN,i) + wr(IDN,i));
    hl = wl(IPR,i)/wl(IDN,i)*(kappal + 1.) + 0.5*(_sqr(wl(ivx,i)) + _sqr(wl(ivy,i)) + _sqr(wl(ivz,i))) + lel;
    hr = wr(IPR,i)/wr(IDN,i)*(kappar + 1.) + 0.5*(_sqr(wr(ivx,i)) + _sqr(wr(ivy,i)) + _sqr(wr(ivz,i))) + ler;
    cbar = sqrt(0.5*gamma*(wl(IPR,i)+ wr(IPR,i))/rhobar);
    pbar = 0.5*(wl(IPR,i) + wr(IPR,i)) + 0.5*(rhobar*cbar)*(wl(ivx,i) - wr(ivx,i));
    ubar = 0.5*(wl(ivx,i) + wr(ivx,i)) + 0.5/(rhobar*cbar)*(wl(IPR,i) - wr(IPR,i));

    if (ubar > 0.) {
      // volume mixing ratio to mass mixing ratio
      Real rd = 1.;
      for (int n = 1; n < ITR; ++n)
        rd -= wl(n,i);

      flx(IDN,i) = ubar*wl(IDN,i)*rd;
      for (int n = 1; n < ITR; ++n)
        flx(n,i) = ubar*wl(IDN,i)*wl(n,i);
      flx(ivx,i) = ubar*wl(IDN,i)*wl(ivx,i) + pbar;
      flx(ivy,i) = ubar*wl(IDN,i)*wl(ivy,i);
      flx(ivz,i) = ubar*wl(IDN,i)*wl(ivz,i);
      flx(IEN,i) = ubar*wl(IDN,i)*hl;
    } else {
      // volume mixing ratio to mass mixing ratio
      Real rd = 1.;
      for (int n = 1; n < ITR; ++n)
        rd -= wr(n,i);

      flx(IDN,i) = ubar*wr(IDN,i)*rd;
      for (int n = 1; n < ITR; ++n)
        flx(n,i) = ubar*wr(IDN,i)*wr(n,i);
      flx(ivx,i) = ubar*wr(IDN,i)*wr(ivx,i) + pbar;
      flx(ivy,i) = ubar*wr(IDN,i)*wr(ivy,i);
      flx(ivz,i) = ubar*wr(IDN,i)*wr(ivz,i);
      flx(IEN,i) = ubar*wr(IDN,i)*hr;
    }
  }
}
