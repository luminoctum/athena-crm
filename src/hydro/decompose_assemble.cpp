// Athena++ headers
#include "hydro.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../eos/eos.hpp"
#include "../coordinates/coordinates.hpp"

void Hydro::DecomposePressure(AthenaArray<Real> &w)
{
  int is = pmy_block->is; int js = pmy_block->js; int ks = pmy_block->ks;
  int ie = pmy_block->ie; int je = pmy_block->je; int ke = pmy_block->ke;
  Coordinates *pco = pmy_block->pcoord;
  EquationOfState *peos = pmy_block->peos;
  Real grav = psrc->GetG1();  // positive in the x-increasing direction
  if (grav == 0.) return;

  // reconstruct top boundary pressure assuming local hydrostatic balance
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      psi_(k,j,ie+1) = w(IPR,k,j,ie) + grav*w(IDN,k,j,ie)*(pco->x1f(ie+1) - pco->x1v(ie));

  // integrate downward from top boundary
  Real gamma = peos->GetGamma();
  Real kappa = 1. - 1./gamma;
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = ie; i >= is; --i) {
        psi_(k,j,i) = psi_(k,j,i+1) - grav*w(IDN,k,j,i)*(pco->x1f(i+1) - pco->x1f(i));
        ps_(k,j,i) = pow(kappa*(psi_(k,j,i+1) - psi_(k,j,i))/(pow(psi_(k,j,i+1),kappa) - pow(psi_(k,j,i),kappa)), gamma);
        w(IPR,k,j,i) -= ps_(k,j,i);
      }

  Real x1min = pmy_block->pmy_mesh->mesh_size.x1min;
  Real x1max = pmy_block->pmy_mesh->mesh_size.x1max;

  // bottom ghost cells
  if (pco->x1v(is-1) < x1min) grav = - psrc->GetG1();
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is-1; i >= 0; --i) {
        psi_(k,j,i) = psi_(k,j,i+1) - grav*w(IDN,k,j,i)*(pco->x1f(i+1) - pco->x1f(i));
        ps_(k,j,i) = pow(kappa*(psi_(k,j,i+1) - psi_(k,j,i))/(pow(psi_(k,j,i+1),kappa) - pow(psi_(k,j,i),kappa)), gamma);
        w(IPR,k,j,i) -= ps_(k,j,i);
      }

  // top ghose cells
  if (pco->x1v(ie+1) > x1max) grav = - psrc->GetG1();
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = ie+1; i < w.GetDim1(); ++i) {
        psi_(k,j,i+1) = psi_(k,j,i) + grav*w(IDN,k,j,i)*(pco->x1f(i+1) - pco->x1f(i));
        ps_(k,j,i) = pow(kappa*(psi_(k,j,i+1) - psi_(k,j,i))/(pow(psi_(k,j,i+1),kappa) - pow(psi_(k,j,i),kappa)), gamma);
        w(IPR,k,j,i) -= ps_(k,j,i);
      }
}

void Hydro::AssemblePressure(AthenaArray<Real> &w)
{
  int js = pmy_block->js; int ks = pmy_block->ks;
  int je = pmy_block->je; int ke = pmy_block->ke;
  Real grav = psrc->GetG1();
  if (grav == 0.) return;

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = 0; i < w.GetDim1(); ++i)
        w(IPR,k,j,i) += ps_(k,j,i);
}
