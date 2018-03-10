#ifndef THERMODYNAMICS_HPP
#define THERMODYNAMICS_HPP
#include "../athena.hpp"

// frequently used inline functions
inline Real SatVaporPresIdeal(Real t, Real p3, Real beta, Real delta) {
  return p3*exp((1. - 1./t)*beta - delta*log(t));
}

// alpha = L/cv evaluated at current temperature
inline Real GasCloudIdeal(Real const prim[], int iv, int ic,
  Real p3, Real t3, Real alpha, Real beta, Real delta)
{
  Real xv = prim[iv];
  Real xc = prim[ic];
  Real t = prim[IDN]/t3;
  Real s = SatVaporPresIdeal(t,p3,beta,delta)/prim[IPR];

  // if saturation vapor pressure is larger than the total pressure
  // evaporate all condensates
  if (s > 1.) return -xc;

  Real g = 1.;
  for (int n = ICD; n < ICD + NVAPOR; ++n) g -= prim[n];
  g -= xv;

  Real s1 = s/(1. - s);
  Real rate = (xv - s1*g)/(1. + alpha*g*(beta/t - delta)*s1/(1. - s));

  // condensate at most xv vapor
  if (rate > 0.) rate = std::min(rate, xv);

  // evaporate at most xc cloud
  if (rate < 0.) rate = std::max(rate, -xc);
  return rate;
}

inline Real f_gas(Real const prim[])
{
  Real fgas = 1.;
  for (int n = ICD; n < ICD + NVAPOR; ++n)
    fgas -= prim[n];
  return fgas;
}

inline Real f_eps(Real const prim[], Real const eps[])
{
  Real feps = 1.;
  for (int n = 1; n < ICD + NVAPOR; ++n)
    feps += prim[n]*(eps[n] - 1.);
  return f_gas(prim)/feps;
}

inline Real f_rcp(Real const prim[], Real const rcp[])
{
  Real frcp= 1.;
  for (int n = 1; n < ICD + NVAPOR; ++n)
    frcp += prim[n]*(rcp[n] - 1.);
  return frcp/f_gas(prim);
}

inline Real dlnTdlnP(Real const prim[], int const isat[], 
  Real const rcp[], Real const beta[], Real const t3[], Real gamma)
{
  Real cp_hat = gamma/(gamma - 1.)*f_rcp(prim, rcp);

  Real xd = 1.;
  for (int n = 1; n < ICD + NVAPOR; ++n)
    xd -= prim[n];

  Real c1 = 0., c2 = 0., c3 = 0.;
  for (int n = 1; n <= NVAPOR; ++n) {
    if (isat[n]) {
      Real latent = beta[n]*t3[n]/prim[IDN] - beta[n+NVAPOR];
      c1 += prim[n]/xd*latent;
      c2 += prim[n]/xd*latent*latent;
    }
    c3 += prim[n]/xd;
  }

  return (1. + c1)/(cp_hat + (c2 + c1*c1)/(1. + c3));
}

void rk4_integrate_z(Real prim[], int isat[], Real const rcp[], Real const eps[],
  Real const beta[], Real const p3[], Real const t3[],
  Real gamma, Real g_ov_Rd, Real dz);

void rk4_integrate_z_adaptive(
  Real prim[], int isat[], Real const rcp[], Real const eps[],
  Real const beta[], Real const p3[], Real const t3[],
  Real gamma, Real g_ov_Rd, Real dz, Real ftol = 1.E-4);

#endif
