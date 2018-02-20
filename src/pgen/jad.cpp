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

enum {iH2O = 1, iNH3 = 2, iH2Oc = 3, iNH3c = 4, iH2Op = 5, iNH3p = 6};
Real grav;

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
  // calculate temperature
  Real p0 = 1.E5;

  for (int i = is; i <= ie; ++i) {
    user_out_var(0,i) = pmicro->Temp(phydro->w, i);
    user_out_var(1,i) = pmicro->Theta(p0, phydro->w, i);
    user_out_var(2,i) = pmicro->Thetav(p0, phydro->w, i);
    user_out_var(3,i) = pmicro->MSE(grav, phydro->w, i);
  }

  /* test whether the result will change if I add back all the condensates
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i)
        for (int n = 1; n <= NVAPOR; ++n) {
          phydro->u(n,k,j,i) += phydro->u(n+NVAPOR,k,j,i);
          phydro->u(n+NVAPOR,k,j,i) = 0.;
        }
  peos->ConservedToPrimitive(phydro->u, phydro->w, pfield->b, phydro->w, pfield->bcc,
    pcoord, is, ie, js, je, ks, ke);*/
}

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  grav = -phydro->psrc->GetG1();
  Real gamma = peos->GetGamma();
  Real Rd = pmicro->GetRd();
  Real cpd = Rd*gamma/(gamma - 1.);

  Real P0 = pin->GetReal("problem", "P0");
  Real T0 = pin->GetReal("problem", "T0");

  Real qH2O = pin->GetReal("problem", "qH2O");
  Real qNH3 = pin->GetReal("problem", "qNH3");

  int i0 = _locate(pcoord->x1v.data(), 0., pcoord->x1v.GetSize());

  std::cout << i0 << std::endl;
  std::cout << pcoord->x1v(i0) << std::endl;

  for (int n = 0; n < NHYDRO; ++n) phydro->w(n,i0) = 0.;
  phydro->w(iH2O,i0) = qH2O;
  phydro->w(iNH3,i0) = qNH3;

  pmicro->DryAdiabat(phydro->w, T0, P0, grav, ks, js, i0, is, ie, false);

  /*Real prim[NHYDRO];
  for (int n = 0; n < NHYDRO; ++n) prim[n] = 0.;
  prim[iNH3] = qNH3;
  prim[iH2O] = qH2O;
  Real cp = pmicro->Cp(prim);

  std::cout << cp << std::endl;

  for (int i = is; i <= ie; ++i) {
    Real x1 = pcoord->x1v(i);
    Real temp = T0 - grav*x1/cp;
    phydro->w(IPR,i) = P0*pow(temp/T0, cp/Rd);
    phydro->w(iH2O,i) = qH2O;
    phydro->w(iH2Oc,i) = 0.0;
    phydro->w(iNH3,i) = qNH3;
    phydro->w(iNH3c,i) = 0.0;

    // calculate virtual temperature
    Real mu = 1.;
    for (int n = 1; n < ITR; ++n)
      mu += phydro->w(n,i)*(1./pmicro->GetMassRatio(n) - 1.);

    phydro->w(IDN,i) = phydro->w(IPR,i)/(Rd*temp*mu);
  }*/

  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je, ks, ke);

  //pmicro->CalculateTemperature(phydro->u);
  //pmicro->SaturationAdjustment(phydro->u);

  /*pmicro->Precipitation(phydro->w, phydro->u, 200., is, ie ,js, je, ks, ke);
  pmicro->Evaporation(phydro->w, phydro->u, 1.);*/


  /* debug
  for (int n = 0; n < NHYDRO; ++n)
    std::cout << phydro->w(n,is) << " ";
  std::cout << std::endl;
  peos->ConservedToPrimitive(phydro->u, phydro->w1, pfield->b,
                                 phydro->w, pfield->bcc, pcoord,
                                 is, ie, js, je, ks, ke);

  for (int n = 0; n < NHYDRO; ++n)
    std::cout << phydro->u(n,(is+ie)/2.) << " ";
  std::cout << std::endl;
  for (int n = 0; n < NHYDRO; ++n)
    std::cout << phydro->w(n,(is+ie)/2.) << " ";
  std::cout << std::endl;*/
}
