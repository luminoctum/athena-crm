//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file wicker.cpp
//  \brief Problem generator for a dry rising bubble
//
// REFERENCE: Wicker and Skamarock, 1998
//========================================================================================

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../math_funcs.hpp"
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../microphysics/microphysics.hpp"
#include "../microphysics/thermodynamics.hpp"

enum {iH2O = 1, iH2Oc = 2};

Real p0, grav;

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  AllocateUserOutputVariables(5);
  SetUserOutputVariableName(0, "temp");
  SetUserOutputVariableName(1, "theta");
  SetUserOutputVariableName(2, "theta_v");
  SetUserOutputVariableName(3, "mse");
  SetUserOutputVariableName(4, "theta_e");
}

void MeshBlock::UserWorkInLoop()
{
  for (int j = js; j <= je; ++j)
    for (int i = is; i <= ie; ++i) {
      user_out_var(0,j,i) = pmicro->Temp(phydro->w, i, j);
      user_out_var(1,j,i) = pmicro->Theta(p0, phydro->w, i, j);
      user_out_var(2,j,i) = pmicro->ThetaV(p0, phydro->w, i, j);
      user_out_var(3,j,i) = pmicro->MSE(grav, phydro->w, i, j);
      user_out_var(4,j,i) = pmicro->ThetaE(p0, phydro->w, i, j);
    }
}

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  grav = - phydro->psrc->GetG1();
  p0 = pin->GetReal("problem", "p0");
  Real Ts = pin->GetReal("problem", "Ts");

  Real xc = pin->GetReal("problem", "xc");
  Real zc = pin->GetReal("problem", "zc");
  Real xr = pin->GetReal("problem", "xr");
  Real zr = pin->GetReal("problem", "zr");
  Real dT = pin->GetReal("problem", "dT");
  Real qt = pin->GetReal("problem", "qt");
  Real p3 = pin->GetReal("microphysics", "p3");
  Real t3 = pin->GetReal("microphysics", "t3");

  Real beta, delta;
  pin->GetReals("microphysics", "beta", &beta, &delta);

  Real Rd = pmicro->GetRd();
  Real gamma = peos->GetGamma();
  Real cpd = gamma/(gamma - 1.)*Rd;

  // moist adiabat
  phydro->w(iH2O,js,is) = qt;
  pmicro->MoistAdiabat(phydro->w, Ts, p0, grav, is, ie, is, js);
  for (int n = 0; n < NHYDRO; ++n)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        phydro->w(n,j,i) = phydro->w(n,js,i);
        phydro->w(IVX,j,i) = 0.;
        phydro->w(IVY,j,i) = 0.;
      }

  // add temperature anomaly
  Real prim[NHYDRO];
  for (int j = js; j <= je; ++j)
    for (int i = is; i <= ie; ++i) {
      Real x1 = pcoord->x1v(i);
      Real x2 = pcoord->x2v(j);
      Real L = sqrt(_sqr((x2 - xc)/xr) + _sqr((x1 - zc)/zr));
      if (L < 1.) {
        pmicro->Hydro2Prim(prim,phydro->w,i,j);
        Real theta = pmicro->Theta(p0, phydro->w,i,j);
        prim[IDN] +=  dT*_sqr(cos(M_PI*L/2.))*pow(prim[IPR]/p0, Rd/cpd)*theta/300.;
        Real rate = GasCloudIdeal(prim, iH2O, iH2Oc, p3, t3, 0., beta, delta);
        prim[iH2O] -= rate;
        prim[iH2Oc] += rate;
        pmicro->Prim2Hydro(prim,phydro->w,i,j);
      }
    }

  UserWorkInLoop();

  /*for (int j = js; j <= je; ++j)
    for (int i = is; i <= ie; ++i) {
      Real x1 = pcoord->x1v(i);
      Real x2 = pcoord->x2v(j);
      Real L = sqrt(_sqr((x2 - xc)/xr) + _sqr((x1 - zc)/zr));
      Real temp = Ts - grav*x1/cpd;
      phydro->w(IPR,j,i) = p0*pow(temp/Ts, cpd/Rd);
      if (L < 1.)
        temp += dT*_sqr(cos(M_PI*L/2.))*pow(phydro->w(IPR,j,i)/p0, Rd/cpd);
      phydro->w(IDN,j,i) = phydro->w(IPR,j,i)/(Rd*temp);
      phydro->w(IVX,j,i) = 0.;
      phydro->w(IVY,j,i) = 0.;
    }*/

  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je, ks, ke);
}
