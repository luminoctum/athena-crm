#include <cstdlib>

#include "../athena_arrays.hpp"
#include "../hydro/hydro.hpp"
#include "../math_funcs.hpp"
#include "microphysics.hpp"

// molecule index
enum {iH2O = 1, iNH3 = 2, iH2Oc = 3, iNH3c = 4, iH2Op = 5, iNH3p = 6};

void Microphysics::SaturationAdjustment(AthenaArray<Real> &u)
{
  std::stringstream msg;
  Real prim[NHYDRO];

  MeshBlock *pmb = pmy_block_;
  Real gamma = pmb->peos->GetGamma();

  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real rho = u(IDN,k,j,i), sum = 1.;
        for (int n = 1; n < ITR; ++n) {
          rho += u(n,k,j,i);
          prim[n] = u(n,k,j,i)/u(IDN,k,j,i)/eps_[n];
          sum += prim[n];
        }
        for (int n = 1; n < ITR; ++n)
          prim[n] /= sum;

        // total energy
        Real u0 = 0., cv = 1.;
        for (int n = 1; n < ICD; ++n)
          u0 -= beta_[n]*t3_[n]*prim[n+NVAPOR];
        for (int n = 1; n < ITR; ++n)
          cv += (rcv_[n]*eps_[n] - 1.)*prim[n];
        u0 += 1./(gamma - 1)*cv*T(k,j,i);

        Real qsig = 1., qeps = 1., qtol = 1.;
        for (int n = 1; n < ICD; ++n) {
          qsig += (prim[n] + prim[n+NVAPOR])*(rcv_[n]*eps_[n] - 1.);
          qeps += (prim[n] + prim[n+NVAPOR])*(eps_[n] - 1.);
          qtol -= prim[n+NVAPOR];
        }

        prim[IDN] = T(k,j,i);
        prim[IPR] = rho*Rd_*T(k,j,i)*qtol/qeps;

        int iter = 0, max_iter = 10;
        // debug
        //std::cout << "Adjustment start, u = " << u(IEN,k,j,i) << std::endl;

        Real Told = 0.;
        while (fabs(prim[IDN] - Told) > 1.E-4 && iter < max_iter) {
          /* debug 
          for (int n = 0; n < ITR; ++n)
            std::cout << prim[n] << " ";
          std::cout << prim[IPR] << " ";
          std::cout <<
          SatVaporPresIdeal(prim[IDN]/t3_[iH2O],p3_[iH2O],beta_[iH2O],beta_[iH2Oc])
          << std::endl;*/

          Real t, alpha, rate, u1;
          
          // condensation
          for (int n = 1; n < 1 + NVAPOR; ++n) {
            int nc = n + NVAPOR;
            t = prim[IDN]/t3_[n];
            alpha = (gamma - 1.)*(beta_[n]/t - beta_[nc] - 1.)/qsig;
            rate = GasCloudIdeal(prim, n, nc,
              p3_[n], t3_[n], alpha, beta_[n], beta_[nc]);
            prim[n] -= rate;
            prim[nc] += rate;
            qtol -= rate;
          }

          // calculate new temperature
          u1 = u0, cv = 1.;
          for (int n = 1; n < ICD; ++n)
            u1 += beta_[n]*t3_[n]*prim[n+NVAPOR];
          for (int n = 1; n < ITR; ++n)
            cv += (rcv_[n]*eps_[n] - 1.)*prim[n];
          Told = prim[IDN];
          prim[IDN] = (gamma - 1.)*u1/cv;

          // adjust pressure
          prim[IPR] = rho*Rd_*prim[IDN]*qtol/qeps;
          iter++;
          if (iter >= max_iter)
            std::cerr << "Saturation Adjustment Iteration reaches maximum."
                      << std::endl;
        }
        // debug
        //std::cout << "Adjustment end, u = ";

        // molar to mass mixing ratios
        sum = 1.;
        for (int n = 1; n < ITR; ++n) {
          sum += prim[n]*(eps_[n] - 1.);
          prim[n] *= eps_[n];
        }
        for (int n = 1; n < ITR; ++n)
          u(n,k,j,i) = rho*prim[n]/sum;

        // update temperature
        T(k,j,i) = prim[IDN];

        /* debug recalculate energy
        Real cvd = Rd_/(gamma - 1);
        Real ue = u(IDN,k,j,i)*cvd*prim[IDN];
        for (int n = 1; n < ITR; ++n)
          ue += u(n,k,j,i)*rcv_[n]*cvd*prim[IDN];
        ue += 0.5/rho*(_sqr(u(IM1,k,j,i)) + _sqr(u(IM2,k,j,i)) + _sqr(u(IM3,k,j,i)));
        for (int n = ICD; n < ICD + NVAPOR; ++n)
          ue += latent_[n]*u(n,k,j,i);
        std::cout << ue << std::endl << std::endl;*/
      }
}

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

  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i)
        for (int n = 0; n < NVAPOR; ++n) {
          Real drho;
          int nc = ICD + n;
          int np = ITR + n;
          if (u(nc,k,j,i) > tiny_number_) {
            if (recondense_(n,k,j,i))
              drho = u(nc,k,j,i);
            else
              drho = u(nc,k,j,i)*dt/autoc_;
            u(nc,k,j,i) -= drho;
            u(np,k,j,i) += drho;
            u(IEN,k,j,i) -= Rd_/(gamma - 1)*drho*rcv_[nc]*T(k,j,i) + latent_[nc]*drho;
          }
        }
}

void Microphysics::Evaporation(AthenaArray<Real> &u, Real dt)
{
  MeshBlock *pmb = pmy_block_;

  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i)
        for (int n = 0; n < NVAPOR; ++n) {
          int ng = 1 + n;
          int nc = ICD + n;
          int np = ITR + n;
          recondense_(n,k,j,i) = false;
          if (u(nc,k,j,i) < tiny_number_ && u(np,k,j,i) > tiny_number_) {
            // limit evaporation rate
            Real h = exp(-t3_[ng]*beta_[ng]/T(k,j,i));
            Real drho = std::min(dt*evapr_*h, u(np,k,j,i));
            u(nc,k,j,i) += drho;
            u(np,k,j,i) -= drho;
            recondense_(n,k,j,i) = true;
          }
        }
}
