#include <cstdlib>

#include "../athena_arrays.hpp"
#include "../hydro/hydro.hpp"
#include "../math_funcs.hpp"
#include "microphysics.hpp"

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
        Real Told = 0.;
        while (fabs(prim[IDN] - Told) > 1.E-4 && iter < max_iter) {
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

        // molar to mass mixing ratios
        sum = 1.;
        for (int n = 1; n < ITR; ++n) {
          sum += prim[n]*(eps_[n] - 1.);
          prim[n] *= eps_[n];
        }
        for (int n = 1; n < ITR; ++n)
          u(n,k,j,i) = rho*prim[n]/sum;

        // update temperature and pressure
        T(k,j,i) = prim[IDN];
        P(k,j,i) = prim[IPR];
      }
}
