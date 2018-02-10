// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../globals.hpp"
#include "../math_funcs.hpp"
#include "../microphysics/microphysics.hpp"
#include "../utils/utils.hpp"

// molecules
enum {iH2O = 1, iNH3 = 2, iH2Oc = 3, iNH3c = 4, iH2Op = 5, iNH3p = 6};
Real friction, heating;

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  AllocateUserOutputVariables(2);
  SetUserOutputVariableName(0, "temp");
  SetUserOutputVariableName(1, "pot_temp");
}

void MeshBlock::UserWorkInLoop()
{
  // calculate temperature
  Real prim[NHYDRO];
  Real p0 = 1.E5;
  Real gamma = peos->GetGamma();

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        for (int n = 0; n < NHYDRO; ++n)
          prim[n] = phydro->w(n,k,j,i);
        Real temp = pmicro->Temperature(prim);
        Real pres = phydro->w(IPR,k,j,i);
        user_out_var(0,k,j,i) = temp;
        user_out_var(1,k,j,i) = temp*pow(p0/pres, (gamma - 1.)/gamma);
      }
}

void JupiterForcing(MeshBlock *pmb, Real const time, Real const dt,
    AthenaArray<Real> const &w, AthenaArray<Real> const &bcc, AthenaArray<Real> &u)
{
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real vx = w(IVX,k,j,i),
             vy = w(IVY,k,j,i),
             vz = w(IVZ,k,j,i);
        Real rho = w(IDN,k,j,i);

        /* turbulent bottom friction
        if (w(IPR,k,j,i) > 150.E5) {
          u(IM1,k,j,i) += - dt*friction*vx*rho; 
          u(IM2,k,j,i) += - dt*friction*vy*rho;
          u(IM3,k,j,i) += - dt*friction*vz*rho;
        }*/
        
        // radiative cooling
        if (w(IPR,k,j,i) < 1.E5) 
          u(IEN,k,j,i) += heating*rho*dt;
      }
}

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  EnrollUserExplicitSourceFunction(JupiterForcing);
}

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  Real grav = -phydro->psrc->GetG1();
  Real gamma = peos->GetGamma();
  Real Rd = pmicro->GetRd();
  Real cpd = Rd*gamma/(gamma - 1.);

  Real P0 = pin->GetReal("problem", "P0");
  Real T0 = pin->GetReal("problem", "T0");

  friction = pin->GetReal("problem", "friction");
  heating = pin->GetReal("problem", "heating");
  Real xH2O = pin->GetReal("problem", "xH2O");
  Real xNH3 = pin->GetReal("problem", "xNH3");

  Real prim[NHYDRO];
  for (int n = 0; n < NHYDRO; ++n) prim[n] = 0.;
  prim[iNH3] = xNH3;
  prim[iH2O] = xH2O;
  Real cp = pmicro->Cp(prim);

  Real kx = 20.*(PI)/(pmy_mesh->mesh_size.x1max - pmy_mesh->mesh_size.x1min);
  Real ky = 20.*(PI)/(pmy_mesh->mesh_size.x2max - pmy_mesh->mesh_size.x2min);
  long int iseed = -1;

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        Real x1 = pcoord->x1v(i);
        Real temp = T0 - grav*x1/cp;
        phydro->w(IPR,k,j,i) = P0*pow(temp/T0, cpd/Rd);
        phydro->w(iH2O,k,j,i) = xH2O;
        phydro->w(iH2Oc,k,j,i) = 0.0;
        phydro->w(iNH3,k,j,i) = xNH3;
        phydro->w(iNH3c,k,j,i) = 0.0;

        // calculate virtual temperature
        Real mu = 1.;
        for (int n = 1; n < ITR; ++n)
          mu += phydro->w(n,k,j,i)*(pmicro->GetMassRatio(n) - 1.);

        phydro->w(IDN,k,j,i) = phydro->w(IPR,k,j,i)*mu/(Rd*temp);
        phydro->w(IVX,k,j,i) = 1.*(ran2(&iseed)-0.5)
          *(1.0+cos(kx*pcoord->x1v(i)))*(1.0+cos(ky*pcoord->x2v(j)))/4.0;
      }

  // saturation adjustment
  pmicro->SaturationAdjustment(phydro->w, is, ie, js, je, ks, ke);

  // remove all condensates
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        phydro->w(iH2Oc,k,j,i) = 0.0;
        phydro->w(iNH3c,k,j,i) = 0.0;
      }

  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je, ks, ke);
}
