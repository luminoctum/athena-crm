#include <cstdlib>

#include "../athena_arrays.hpp"
#include "../hydro/hydro.hpp"
#include "../math_funcs.hpp"
#include "microphysics.hpp"
#include "thermodynamics.hpp"

void Microphysics::Evaporation(AthenaArray<Real> &u, Real dt)
{
  MeshBlock *pmb = pmy_block_;
  Real gamma_d = pmb->peos->GetGamma();
  Real tiny_number = 1.0E-20;

  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        // calculate total density
        Real rho = 0.;
        for (int n = 0; n < ITR; ++n)
          rho += u(n,k,j,i);

        for (int n = 0; n < NVAPOR; ++n) {
          int ng = 1 + n;
          int nc = ICD + n;
          int np = ITR + n;
          recondense_(n,k,j,i) = false;

          // This subroutine is the first microphysics step executed after
          // the dynamics step. It is possible that some variables might be negative
          // after dynamics.
          // set negative value to zero
          u(ng,k,j,i) = std::max(0., u(ng,k,j,i));
          u(nc,k,j,i) = std::max(0., u(nc,k,j,i));
          u(np,k,j,i) = std::max(0., u(np,k,j,i));
          if (u(nc,k,j,i) > tiny_number) continue;  // don't evaporate when there is cloud

          Real feps = 1.;
          for (int n = ICD; n < ICD + NVAPOR; ++n)
            feps -= u(n,k,j,i)/rho;
          for (int n = 1; n < 1 + NVAPOR; ++n)
            feps += u(n,k,j,i)/rho*(1./eps_[n] - 1.);

          Real gamma = rcp_[ng]/rcv_[ng]*gamma_d;
          Real t = temp_(k,j,i)/t3_[ng];
          Real svp = SatVaporPresIdeal(t,p3_[ng],beta_[ng],beta_[nc]);
          Real qs = svp*eps_[ng]/pres_(k,j,i)*feps;
          Real qa = u(ng,k,j,i)/rho;
          if (qa > qs) continue;  // saturated air, usually this is not used

          Real rate, LovRT, dmdt;
          if (qs > 1.) // boiling
            rate = u(np,k,j,i)/dt;
          else {
            LovRT = beta_[ng]/t - beta_[nc];
            dmdt = Kw_/_sqr(Dp_)*(qs - qa)/(1. + Kw_/Kt_*qs*_sqr(LovRT)/(beta_[nc] + gamma/(gamma - 1)));
            rate = 12.*u(np,k,j,i)*rho/rhol_*dmdt;
          }
          //std::cout << qs << " " << qa << " " << LovRT << " " << dmdt << std::endl;
          Real drho = std::min(dt*rate, u(np,k,j,i));
          u(nc,k,j,i) += drho;
          u(np,k,j,i) -= drho;
          u(IEN,k,j,i) += Rd_/(gamma_d - 1)*drho*rcv_[nc]*temp_(k,j,i) - latent_[nc]*drho;
          recondense_(n,k,j,i) = true;
        }
      }
}
