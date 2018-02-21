#include <cstdlib>

#include "../athena_arrays.hpp"
#include "../hydro/hydro.hpp"
#include "../math_funcs.hpp"
#include "microphysics.hpp"

void Microphysics::Evaporation(AthenaArray<Real> &u, Real dt)
{
  MeshBlock *pmb = pmy_block_;
  Real gamma = pmb->peos->GetGamma();

  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        // calculate partial pressure of dry air
        Real sum = 1.;
        for (int n = 1; n < 1 + NVAPOR; ++n)
          sum += u(n,k,j,i)/u(IDN,k,j,i)/eps_[n];
        Real pd = 1./sum*P(k,j,i);

        for (int n = 0; n < NVAPOR; ++n) {
          int ng = 1 + n;
          int nc = ICD + n;
          int np = ITR + n;
          recondense_(n,k,j,i) = false;

          Real t = T(k,j,i)/t3_[ng];
          Real qs = SatVaporPresIdeal(t, p3_[ng], beta_[ng], beta_[nc])/pd*eps_[ng];
          Real qa = u(ng,k,j,i)/u(IDN,k,j,i);

          Real rate = (qs - qa)/(1. + (gamma - 1)/gamma*_sqr(beta_[ng]/t - beta_[nc])*qs);
          Real drho = std::min(dt*evapr_*rate, u(np,k,j,i));
          //Real drho = u(np,k,j,i);
          u(nc,k,j,i) += drho;
          u(np,k,j,i) -= drho;
          u(IEN,k,j,i) += Rd_/(gamma - 1)*drho*rcv_[nc]*T(k,j,i) + latent_[nc]*drho;
          recondense_(n,k,j,i) = true;

        }
      }
}
