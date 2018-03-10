#include <iostream>
#include <cstdlib>

#include "../athena_arrays.hpp"
#include "../hydro/hydro.hpp"
#include "../math_funcs.hpp"
#include "microphysics.hpp"
#include "thermodynamics.hpp"

void Microphysics::MoistAdiabat(AthenaArray<Real>& w, Real T0, Real P0, Real grav,
  int k, int j, int i0, int is, int ie, Real ptop, Real pbot) const
{
  Coordinates *pcoord = pmy_block_->pcoord;
  Real gamma = pmy_block_->peos->GetGamma();

  // mass to molar mixing ratio
  Real prim1[NHYDRO], prim2[NHYDRO], sum = 1.;
  for (int n = 1; n < ITR; ++n) {
    prim1[n] = w(n,k,j,i0)/eps_[n];
    sum += w(n,k,j,i0)*(1./eps_[n] - 1.);
  }
  for (int n = 1; n < ITR; ++n)
    prim1[n] /= sum;
  prim1[IDN] = T0;
  prim1[IPR] = P0;

  for (int n = 0; n < NHYDRO; ++n)
    prim2[n] = prim1[n];

  Real rcp[ITR];
  for (int n = 0; n < ITR;  ++n)
    rcp[n] = rcp_[n]*eps_[n];

  int isat1[1+NVAPOR], isat2[1+NVAPOR];
  for (int n = 0; n <= NVAPOR; ++n)
    isat1[n] = 0;

  // equilibrate at i0
  for (int n = 1; n < 1 + NVAPOR; ++n) {
    int nc = n + NVAPOR;
    Real rate = GasCloudIdeal(prim1, n, nc, p3_[n], t3_[n], 0., beta_[n], beta_[nc]);
    //Real rate = 0.;
    prim1[n] -= rate;
    prim1[nc] += rate;
    if (rate > 0.) isat1[n] = 1;
  }
  for (int n = 0; n <= NVAPOR; ++n)
    isat2[n] = isat1[n];
  SetPrimitive(prim1, w, i0, j, k);

  // upward
  for (int i = i0+1; i <= ie; ++i) {
    // RK4 integration 
    rk4_integrate_z_adaptive(prim1, isat1, rcp, eps_, beta_, p3_, t3_, gamma,
      grav/Rd_, pcoord->dx1v(i));
    SetPrimitive(prim1, w, i, j, k);
  }

  // downward
  for (int i = i0-1; i >= is; --i) {
    // RK4 integration 
    rk4_integrate_z_adaptive(prim2, isat2, rcp, eps_, beta_, p3_, t3_, gamma,
      grav/Rd_, -pcoord->dx1v(i));
    SetPrimitive(prim2, w, i, j, k);
  }
}
