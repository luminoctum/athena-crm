#include <cstdlib>

#include "../athena_arrays.hpp"
#include "../hydro/hydro.hpp"
#include "../math_funcs.hpp"
#include "microphysics.hpp"
#include "thermodynamics.hpp"

void Microphysics::DryAdiabat(AthenaArray<Real>& w, Real T0, Real P0, Real grav,
  int k, int j, int i0, int is, int ie, Real ptop, Real pbot) const
{
  Coordinates *pcoord = pmy_block_->pcoord;
  // mass to molar mixing ratio
  Real prim[NHYDRO], prim0[NHYDRO], sum = 1.;
  for (int n = 1; n < ITR; ++n) {
    prim[n] = w(n,k,j,i0)/eps_[n];
    sum += w(n,k,j,i0)*(1./eps_[n] - 1.);
  }
  for (int n = 1; n < ITR; ++n)
    prim[n] /= sum;
  prim[IDN] = T0;
  prim[IPR] = P0;

  for (int n = 0; n < NHYDRO; ++n)
    prim0[n] = prim[n];

  // upward
  Real chi, Tv, Tv0, rate;
  for (int i = i0; i <= ie; ++i) {
    // condensation at this pressure
    // debug
    //for (int n = 0; n < NHYDRO; ++n)
    //  std::cout << prim[n] << " ";
    //std::cout << std::endl;

    for (int n = 1; n < 1 + NVAPOR; ++n) {
      int nc = n + NVAPOR;
      rate = GasCloudIdeal(prim, n, nc, p3_[n], t3_[n], 0., beta_[n], beta_[nc]);
      prim[n] -= rate;
      prim[nc] += rate;
    }

    // debug
    //for (int n = 0; n < NHYDRO; ++n)
    //  std::cout << prim[n] << " ";
    //std::cout << std::endl;

    // remove condensates
    for (int n = ICD; n < ICD + NVAPOR; ++n)
      prim[n] = 0.;

    // set mass mixing ratio
    sum = 1.;
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

    Tv0 = prim[IDN]*feps;
    w(IPR,k,j,i) = prim[IPR];
    w(IDN,k,j,i) = prim[IPR]/(Rd_*Tv0);

    chi = Chi(w,i,j,k);
    if (prim[IPR] < ptop) {
      Tv = Tv0;
      prim[IPR] = w(IPR,k,j,i)*exp(-grav*pcoord->dx1v(i)/(Rd_*Tv));
    } else {
      Tv = Tv0 - chi*grav*pcoord->dx1v(i)/Rd_;
      prim[IPR] = w(IPR,k,j,i)*pow(Tv/Tv0, 1./chi);
    }

    // set tp for next grid
    prim[IDN] = Tv/feps;
  }

  // downward
  for (int i = i0; i > is; --i) {
    chi = Chi(w,i,j,k);
    Tv0 = Tempv(w,i,j,k);
    if (w(IPR,k,j,i) > pbot) {
      Tv = Tv0;
      prim[IPR] = w(IPR,k,j,i)*exp(grav*pcoord->dx1v(i)/(Rd_*Tv));
    } else {
      Tv = Tv0 + chi*grav*pcoord->dx1v(i)/Rd_;
      prim[IPR] = w(IPR,k,j,i)*pow(Tv/Tv0, 1./chi);
    }
    prim[IDN] = Tv/Tv0*Temp(w,i,j,k);

    // set initial mixing ratio
    for (int n = 1; n < ITR; ++n)
      prim[n] = prim0[n];

    // condensation at this pressure
    for (int n = 1; n < 1 + NVAPOR; ++n) {
      int nc = n + NVAPOR;
      rate = GasCloudIdeal(prim, n, nc, p3_[n], t3_[n], 0., beta_[n], beta_[nc]);
      prim[n] -= rate;
      prim[nc] += rate;
    }

    // remove condensates
    for (int n = ICD; n < ICD + NVAPOR; ++n)
      prim[n] = 0.;

    // set mass mixing ratio
    sum = 1.;
    for (int n = 1; n < ITR; ++n) {
      w(n,k,j,i-1) = prim[n]*eps_[n];
      sum += prim[n]*(eps_[n] - 1.); 
    }
    for (int n = 1; n < ITR; ++n)
      w(n,k,j,i-1) /= sum;

    // virtual temperature
    Real feps = 1.;
    for (int n = ICD; n < ICD + NVAPOR; ++n)
      feps -= w(n,k,j,i-1);
    for (int n = 1; n < 1 + NVAPOR; ++n)
      feps += w(n,k,j,i-1)*(1./eps_[n] - 1.);

    w(IPR,k,j,i-1) = prim[IPR];
    w(IDN,k,j,i-1) = prim[IPR]/(Rd_*prim[IDN]*feps);
  }
}
