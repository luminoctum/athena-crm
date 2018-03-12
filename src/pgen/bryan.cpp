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

enum {iH2O = 1, iH2Oc = 2};

Real p0, grav;

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  AllocateUserOutputVariables(5);
  SetUserOutputVariableName(0, "T");
  SetUserOutputVariableName(1, "Theta");
  SetUserOutputVariableName(2, "ThetaV");
  SetUserOutputVariableName(3, "MSE");
  SetUserOutputVariableName(4, "ThetaE");
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

  //Real gamma = peos->GetGamma();
  //Real cpd = gamma/(gamma - 1.)*Rd;

  phydro->w(iH2O,ks,js,is) = qt;
  pmicro->MoistAdiabat(phydro->w, Ts, p0, grav, ks, js, is, is, ie);

  for (int j = js; j <= je; ++j)
    for (int i = is; i <= ie; ++i) {
      Real x1 = pcoord->x1v(i);
      Real x2 = pcoord->x2v(j);
      Real L = sqrt(_sqr((x2 - xc)/xr) + _sqr((x1 - zc)/zr));
      for (int n = 0; n < NHYDRO; ++n)
        phydro->w(n,j,i) = phydro->w(n,js,i);
      //Real temp = pmicro->Temp(phydro->w,i,j,k);
      //if (L < 1.)
      //  temp += dT*_sqr(cos(M_PI*L/2.))*pow(phydro->w(IPR,j,i)/p0, Rd/cpd);
      //phydro->w(IDN,j,i) = phydro->w(IPR,j,i)/(Rd*temp);
      //phydro->w(IVX,j,i) = 0.;
      //phydro->w(IVY,j,i) = 0.;
    }

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

  UserWorkInLoop();
  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je, ks, ke);
}
