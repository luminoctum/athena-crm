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
Real friction, heating, grav, P0, T0, ptop, pbot, diag_dt, next_diag_time;

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
  Real time = pmy_mesh->time;

  // only calculate diagnostic variable at output times
  if ((time == pmy_mesh->start_time) ||
     (time >= pmy_mesh->tlim) ||
     (time >= next_diag_time)) {
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j)
        for (int i = is; i <= ie; ++i) {
          user_out_var(0,k,j,i) = pmicro->Temp(phydro->w,i,j,k);
          user_out_var(1,k,j,i) = pmicro->Theta(P0,phydro->w,i,j,k);
          user_out_var(2,k,j,i) = pmicro->ThetaV(P0,phydro->w,i,j,k);
          user_out_var(3,k,j,i) = pmicro->MSE(grav,phydro->w,i,j,k);
        }
    next_diag_time += diag_dt;
  }
}

void JupiterForcing(MeshBlock *pmb, Real const time, Real const dt,
    AthenaArray<Real> const &w, AthenaArray<Real> const &bcc, AthenaArray<Real> &u)
{
  Real dp = pmb->pmicro->GetDp();

  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j) {
      //u(IM1,k,j,i) += - dt*friction*w(IV1,k,j,i)*w(IDN,k,j,i); 
      //u(IM2,k,j,pmb->is) -= dt*w(IV2,k,j,pmb->is)*w(IDN,k,j,pmb->is)/friction;

      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real v1 = w(IV1,k,j,i),
             v2 = w(IV2,k,j,i),
             v3 = w(IV3,k,j,i);
        Real rho = w(IDN,k,j,i);
        //Real tv = pmb->pmicro->TerminalVelocity(dp, rho, grav);
        Real tv = -20.;

        // radiative cooling
        if (w(IPR,k,j,i) < 0.5E5 && w(IPR,k,j,i) > ptop) 
          u(IEN,k,j,i) += heating*rho*dt;

        // falling precipitations
        Real qp = 0.;
        for (int n = 0; n < NVAPOR; ++n)
          qp += w(ITR + n,k,j,i);
        u(IM1,k,j,i) += -dt*rho*qp*grav;
        u(IEN,k,j,i) += -dt*rho*qp*grav*tv;

        // sponge layerprim
        //u(IM1,k,j,i) -= dt*rho*v1*pow(p0/w(IPR,k,j,i),2)/friction;
        //u(IM2,k,j,i) -= dt*rho*v2*pow(p0/w(IPR,k,j,i),2)/friction;
      }
    }
}

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  grav = - pin->GetReal("hydro", "grav_acc1");
  friction = pin->GetReal("problem", "friction");
  heating = pin->GetReal("problem", "heating");
  P0 = pin->GetReal("problem", "P0");
  T0 = pin->GetReal("problem", "T0");
  ptop = pin->GetReal("problem", "ptop");
  pbot = pin->GetReal("problem", "pbot");
  diag_dt = FLT_MAX;

  // get diagnostic output frequency
  InputBlock *pib = pin->pfirst_block;
  while (pib != NULL) {
    if (pib->block_name.compare(0,6,"output") == 0) {
      try {
        if (pin->GetString(pib->block_name, "variable") == "uov") {
          diag_dt = pin->GetReal(pib->block_name, "dt");
          try {
            next_diag_time = pin->GetReal(pib->block_name, "next_time");
          } catch (std::exception const& ex) {
            next_diag_time = time + diag_dt;
          }
        }
      } catch (std::exception const& ex) {}
    }
    pib = pib->pnext;
  }

  EnrollUserExplicitSourceFunction(JupiterForcing);
}

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  Real qH2O = pin->GetReal("problem", "qH2O");
  Real qNH3 = pin->GetReal("problem", "qNH3");

  int i0 = _locate(pcoord->x1v.data(), 0., pcoord->x1v.GetSize());

  for (int n = 0; n < NHYDRO; ++n)
    phydro->w(n,ks,js,i0) = 0.;
  phydro->w(iH2O,ks,js,i0) = qH2O;
  phydro->w(iNH3,ks,js,i0) = qNH3;

  //pmicro->DryAdiabat(phydro->w, T0, P0, grav, ks, js, i0, is, ie, ptop, pbot);
  pmicro->DryAdiabat(phydro->w, T0, P0, grav, is, ie, i0, js, ks, ptop);

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
      }

  UserWorkInLoop();
  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord,
    is, ie, js, je, ks, ke);
}
