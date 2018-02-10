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
  Real grav = -phydro->psrc->GetG1();
  Real gamma = peos->GetGamma();
  Real Rd = pmicro->GetRd();
  Real cpd = Rd*gamma/(gamma - 1.);

  Real P0 = pin->GetReal("problem", "P0");
  Real T0 = pin->GetReal("problem", "T0");

  Real xH2O = pin->GetReal("problem", "xH2O");
  Real xNH3 = pin->GetReal("problem", "xNH3");

  Real prim[NHYDRO];
  for (int n = 0; n < NHYDRO; ++n) prim[n] = 0.;
  prim[iNH3] = xNH3;
  prim[iH2O] = xH2O;
  Real cp = pmicro->Cp(prim);

  for (int i = is; i <= ie; ++i) {
    Real x1 = pcoord->x1v(i);
    Real temp = T0 - grav*x1/cp;
    phydro->w(IPR,i) = P0*pow(temp/T0, cp/Rd);
    phydro->w(iH2O,i) = xH2O;
    phydro->w(iH2Oc,i) = 0.0;
    phydro->w(iNH3,i) = xNH3;
    phydro->w(iNH3c,i) = 0.0;

    // calculate virtual temperature
    Real mu = 1.;
    for (int n = 1; n < ITR; ++n)
      mu += phydro->w(n,i)*(pmicro->GetMassRatio(n) - 1.);

    phydro->w(IDN,i) = phydro->w(IPR,i)*mu/(Rd*temp);
  }

  phydro->w(iH2Op,is) = 0.0;

  pmicro->SaturationAdjustment(phydro->w, is, ie, js, je, ks, ke);

  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je, ks, ke);

  pmicro->Precipitation(phydro->w, phydro->u, 200., is, ie ,js, je, ks, ke);
  pmicro->Evaporation(phydro->w, phydro->u, 1.);


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
  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je, ks, ke);
}
