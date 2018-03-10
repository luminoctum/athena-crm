#include <iostream>
#include <algorithm>
#include "thermodynamics.hpp"

void rk4_integrate_z(Real prim[], int isat[], Real const rcp[], Real const eps[],
  Real const beta[], Real const p3[], Real const t3[],
  Real gamma, Real g_ov_Rd, Real dz)
{
  Real step[] = {0.5, 0.5, 1.};
  Real temp = prim[IDN];
  Real pres = prim[IPR];
  Real dTdz[4], chi[4];

  for (int rk = 0; rk < 4; ++rk) {
    // reset vapor
    for (int n = 1; n < 1 + NVAPOR; ++n) {
      prim[n] += prim[n+NVAPOR];
      prim[n+NVAPOR] = 0.;
      isat[n] = 0;
    }

    for (int n = 1; n < 1 + NVAPOR; ++n) {
      int nc = n + NVAPOR;
      Real rate = GasCloudIdeal(prim, n, nc, p3[n], t3[n], 0., beta[n], beta[nc]);
      //Real rate = 0.;
      prim[n] -= rate;
      prim[nc] += rate;
      if (rate > 0.) isat[n] = 1;
    }

    // calculate tendency
    chi[rk] = dlnTdlnP(prim, isat, rcp, beta, t3, gamma);
    dTdz[rk] = - chi[rk]*g_ov_Rd/f_eps(prim, eps);
    
    // integrate over dz
    if (rk < 3) {
      prim[IDN] = temp + dTdz[rk]*dz*step[rk];
      prim[IPR] = pres*pow(prim[IDN]/temp, 1./chi[rk]);
    } else {
      prim[IDN] = temp + 1./6.*(dTdz[0] + 2.*dTdz[1] + 2.*dTdz[2] + dTdz[3])*dz;
      Real chi_avg = 1./6.*(chi[0] + 2.*chi[1] + 2.*chi[2] + chi[3]);
      prim[IPR] = pres*pow(prim[IDN]/temp, 1./chi_avg);
    }
  }

  // recondensation
  for (int n = 1; n < 1 + NVAPOR; ++n) {
    int nc = n + NVAPOR;
    Real rate = GasCloudIdeal(prim, n, nc, p3[n], t3[n], 0., beta[n], beta[nc]);
    //Real rate = 0.;
    prim[n] -= rate;
    prim[nc] += rate;
    if (rate > 0.) isat[n] = 1;
  }
}

void rk4_integrate_z_adaptive(Real prim[], int isat[], Real const rcp[], Real const eps[],
  Real const beta[], Real const p3[], Real const t3[],
  Real gamma, Real g_ov_Rd, Real dz, Real ftol)
{
  Real prim1[NHYDRO], prim2[NHYDRO];
  int  isat1[1+NVAPOR], isat2[1+NVAPOR];

  for (int n = 0; n < NHYDRO; ++n) {
    prim1[n] = prim[n];
    prim2[n] = prim[n];
  }
  for (int n = 0; n < 1+ NVAPOR; ++n) {
    isat1[n] = isat[n];
    isat2[n] = isat[n];
  }

  // trail step
  rk4_integrate_z(prim1, isat1, rcp, eps, beta, p3, t3, gamma, g_ov_Rd, dz);

  // refined step
  rk4_integrate_z(prim2, isat2, rcp, eps, beta, p3, t3, gamma, g_ov_Rd, dz/2.);
  rk4_integrate_z(prim2, isat2, rcp, eps, beta, p3, t3, gamma, g_ov_Rd, dz/2.);

  if (fabs(prim2[IDN] - prim1[IDN]) > ftol) {
    rk4_integrate_z_adaptive(prim, isat, rcp, eps, beta, p3, t3, gamma, g_ov_Rd, dz/2., ftol/2.);
    rk4_integrate_z_adaptive(prim, isat, rcp, eps, beta, p3, t3, gamma, g_ov_Rd, dz/2., ftol/2.);
  } else {
    for (int n = 0; n < NHYDRO; ++n)
      prim[n] = prim2[n];
    for (int n = 0; n < 1+ NVAPOR; ++n)
      isat[n] = isat2[n];
  }
}
