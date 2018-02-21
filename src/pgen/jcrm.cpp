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
Real friction, heating, grav, P0;

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  AllocateUserOutputVariables(4);
  SetUserOutputVariableName(0, "temp");
  SetUserOutputVariableName(1, "theta");
  SetUserOutputVariableName(2, "thetav");
  SetUserOutputVariableName(3, "MSE");
}

void MeshBlock::UserWorkInLoop()
{
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        user_out_var(0,k,j,i) = pmicro->Temp(phydro->w,i,j,k);
        user_out_var(1,k,j,i) = pmicro->Theta(P0,phydro->w,i,j,k);
        user_out_var(2,k,j,i) = pmicro->Thetav(P0,phydro->w,i,j,k);
        user_out_var(3,k,j,i) = pmicro->MSE(grav,phydro->w,i,j,k);
      }
}

void JupiterForcing(MeshBlock *pmb, Real const time, Real const dt,
    AthenaArray<Real> const &w, AthenaArray<Real> const &bcc, AthenaArray<Real> &u)
{
  Real p0 = 0.5E5;
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j) {
      //u(IM1,k,j,i) += - dt*friction*w(IV1,k,j,i)*w(IDN,k,j,i); 
      //u(IM2,k,j,pmb->is) -= dt*w(IV2,k,j,pmb->is)*w(IDN,k,j,pmb->is)/friction;

      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real v1 = w(IV1,k,j,i),
             v2 = w(IV2,k,j,i),
             v3 = w(IV3,k,j,i);
        Real rho = w(IDN,k,j,i);

        // radiative cooling
        if (w(IPR,k,j,i) < p0) 
          u(IEN,k,j,i) += heating*rho*dt;

        // sponge layer
        //u(IM1,k,j,i) -= dt*rho*v1*pow(p0/w(IPR,k,j,i),2)/friction;
        //u(IM2,k,j,i) -= dt*rho*v2*pow(p0/w(IPR,k,j,i),2)/friction;
      }
    }
}

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  EnrollUserExplicitSourceFunction(JupiterForcing);
}

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  grav = -phydro->psrc->GetG1();
  friction = pin->GetReal("problem", "friction");
  heating = pin->GetReal("problem", "heating");
  P0 = pin->GetReal("problem", "P0");

  Real T0 = pin->GetReal("problem", "T0");
  Real qH2O = pin->GetReal("problem", "qH2O");
  Real qNH3 = pin->GetReal("problem", "qNH3");

  int i0 = _locate(pcoord->x1v.data(), 0., pcoord->x1v.GetSize());

  for (int n = 0; n < NHYDRO; ++n)
    phydro->w(n,ks,js,i0) = 0.;
  phydro->w(iH2O,ks,js,i0) = qH2O;
  phydro->w(iNH3,ks,js,i0) = qNH3;

  pmicro->DryAdiabat(phydro->w, T0, P0, grav, ks, js, i0, is, ie, false);

  Real kx = 20.*(PI)/(pmy_mesh->mesh_size.x1max - pmy_mesh->mesh_size.x1min);
  Real ky = 20.*(PI)/(pmy_mesh->mesh_size.x2max - pmy_mesh->mesh_size.x2min);
  long int iseed = -1;

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        for (int n = 0; n < NHYDRO; ++n)
          phydro->w(n,k,j,i) = phydro->w(n,ks,js,i);

        // add pertubation
        phydro->w(IVX,k,j,i) = 1.*(ran2(&iseed)-0.5)
          *(1.0+cos(kx*pcoord->x1v(i)))*(1.0+cos(ky*pcoord->x2v(j)))/4.0;

        // diagnostic variables
        user_out_var(0,k,j,i) = pmicro->Temp(phydro->w,i,j,k);
        //std::cout << user_out_var(0,k,j,i) << std::endl;
        user_out_var(1,k,j,i) = pmicro->Theta(P0,phydro->w,i,j,k);
        user_out_var(2,k,j,i) = pmicro->Thetav(P0,phydro->w,i,j,k);
        user_out_var(3,k,j,i) = pmicro->MSE(grav,phydro->w,i,j,k);
      }

  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord,
    is, ie, js, je, ks, ke);
}
