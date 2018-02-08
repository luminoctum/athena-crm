#include <cstdlib>

#include "../athena_arrays.hpp"
#include "../hydro/hydro.hpp"
#include "../math_funcs.hpp"
#include "microphysics.hpp"
#include "../utils/logging_info.hpp"  // for debug

inline Real SatVaporPresIdeal(Real t, Real p3, Real beta, Real gamma) {
  return p3*exp((1. - 1./t)*beta - gamma*log(t));
}

inline Real SatVaporPresH2OIdeal(Real T)
{
  Real betal = 24.88,
       gammal = 5.06,
       betas = 22.98,
       gammas = 0.52,
       t3 = 273.16,
       p3 = 611.7;

  //return T > t3 ? SatVaporPresIdeal(T/t3, p3, betal, gammal)
  return SatVaporPresIdeal(T/t3, p3, betas, gammas);
}

inline Real SatVaporPresNH3Ideal(Real T)
{
  Real betal = 20.08,
       gammal = 5.62,
       betas = 20.64,
       gammas = 1.43,
       t3 = 195.4,
       p3 = 6060.;

  //return T > t3 ? SatVaporPresIdeal(T/t3, p3, betal, gammal)
  return SatVaporPresIdeal(T/t3, p3, betas, gammas);
}

// alpha = L/cv evaluated at current temperature
inline Real GasCloudIdeal(Real const prim[], int iv, int ic, Real temp,
  Real p3, Real t3, Real alpha, Real beta, Real gamma)
{
  Real xv = prim[iv];
  Real xc = prim[ic];
  Real s = SatVaporPresIdeal(temp/t3,p3,beta,gamma)/prim[IPR];

  // if saturation vapor pressure is larger than the total pressure
  // evaporate all condensates
  if (s > 1.) return -xc;

  Real xt = 1.;
  for (int n = ICD; n < ICD + NVAPOR; ++n) xt -= prim[n];

  Real dsdt = s/temp*(beta/temp*t3 - gamma);
  Real rate = (xv - s*xt)/((1. - s)*(1. + alpha*(xt - xv)/_sqr(1. - s)*dsdt));

  // condensate at most xv vapor
  if (rate > 0.) rate = std::min(rate, xv);

  // evaporate at most xc cloud
  if (rate < 0.) rate = std::max(rate, -xc);
  return rate;
}

void Microphysics::SaturationAdjust(AthenaArray<Real> &w, AthenaArray<Real> &u)
{
  std::stringstream msg;
  Real prim[NHYDRO];

  int iH2O = 1;
  int iNH3 = 2;
  int iH2Oc = 3;
  int iNH3c = 4;

  Real p3_H2O = 611.7;
  Real t3_H2O = 273.16;
  Real beta_H2O = 22.98;
  Real gamma_H2O = 0.52;

  Real p3_NH3 = 6060.;
  Real t3_NH3 = 195.4;
  Real beta_NH3 = 20.64;
  Real gamma_NH3 = 1.43;

  MeshBlock *pmb = pmy_block_;
  Real gamma = pmb->peos->GetGamma();

  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        for (int n = 0; n < NHYDRO; ++n) prim[n] = w(n,k,j,i);

        /*std::cout << "before adjust ";
        for (int n = 0; n < NHYDRO; ++n)
          std::cout << prim[n] << " ";
        std::cout << std::endl;*/

        Real gmix = 1., qmass = 1.;
        for (int n = 1; n < ITR; ++n) {
          gmix += prim[n]*(rcv_[n] - 1.);
          qmass += prim[n]*(eps_[n] - 1.);
        }
        Real alpha_H2O = -(gamma - 1.)*latent_[iH2O]*eps_[iH2O]/(Rd_*gmix);
        Real alpha_NH3 = -(gamma - 1.)*latent_[iNH3]*eps_[iNH3]/(Rd_*gmix);

        Real dT = 1., qtol = 1., rate;
        for (int n = ICD; n < ICD + NVAPOR; ++n)
          qtol -= prim[n];

        int iter = 0, max_iter = 10;
        LoggingInfo<double> log1("dT");
        while (fabs(dT) > 1.E-4) {
          Real temp = Temperature(prim);
          // H2O - H2O(s)
          rate = GasCloudIdeal(prim, iH2O, iH2Oc, temp,
            p3_H2O, t3_H2O, alpha_H2O, beta_H2O, gamma_H2O);
          prim[iH2O] -= rate;
          prim[iH2Oc] += rate;
          dT = alpha_H2O*rate;
          qtol -= rate;

          // NH3 - NH3(s)
          rate = GasCloudIdeal(prim, iNH3, iNH3c, temp + dT,
            p3_NH3, t3_NH3, alpha_NH3, beta_NH3, gamma_NH3);
          prim[iNH3] -= rate;
          prim[iNH3c] += rate;
          dT += alpha_NH3*rate;
          qtol -= rate;

          // adjust pressure
          prim[IPR] = prim[IDN]*Rd_*(temp + dT)*qtol/qmass;
          iter++;
          if (iter > max_iter) {
            msg << "Saturation Adjustment Iteration reaches maximum." << std::endl;
            throw std::runtime_error(msg.str().c_str());
          }
          log1 << temp << dT << SatVaporPresH2OIdeal(temp+dT);
        }

        Real qd = 1.;
        for (int n = 1; n < ITR; ++n) qd -= prim[n];

        // adjust conserved variables
        u(iH2O,k,j,i) = u(IDN,k,j,i)*prim[iH2O]/qd*eps_[iH2O];
        u(iH2Oc,k,j,i) = u(IDN,k,j,i)*prim[iH2Oc]/qd*eps_[iH2O];
        u(iNH3,k,j,i) = u(IDN,k,j,i)*prim[iNH3]/qd*eps_[iNH3];
        u(iNH3c,k,j,i) = u(IDN,k,j,i)*prim[iNH3c]/qd*eps_[iNH3];

        // adjust primitive variables
        for (int n = 0; n < NHYDRO; ++n) w(n,k,j,i) = prim[n];

        /*std::cout << "after adjust ";
        for (int n = 0; n < NHYDRO; ++n)
          std::cout << prim[n] << " ";
        std::cout << std::endl;

        if (iter > 1) {
          std::cout << "Number of iterations: " << iter << std::endl;
          std::cout << "dT = " << log1 << std::endl;
        }*/
      }
}
