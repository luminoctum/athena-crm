//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file robert.cpp
//  \brief Problem generator for rising bubble
//
// REFERENCE: Robert et al., 1992
//========================================================================================

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../microphysics/microphysics.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../utils/utils.hpp"
#include "../globals.hpp"

Real grav, p0, Ts, Rd, cp, K;

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  grav = - phydro->psrc->GetG1();
  Rd = pmicro->GetRd();
  p0 = pin->GetReal("problem", "p0");
  Ts = pin->GetReal("problem", "Ts");
  cp = pin->GetReal("problem", "cp");

  Real xc = pin->GetReal("problem", "xc");
  Real zc = pin->GetReal("problem", "zc");
  Real s = pin->GetReal("problem", "s");
  Real a = pin->GetReal("problem", "a");
  Real dT = pin->GetReal("problem", "dT");

  Real qball = pin->GetReal("problem", "qball");

  for (int j = js; j <= je; ++j)
    for (int i = is; i <= ie; ++i) {
      Real x1 = pcoord->x1v(i);
      Real x2 = pcoord->x2v(j);
      Real r = sqrt((x2-xc)*(x2-xc) + (x1-zc)*(x1-zc));
      Real temp = Ts - grav*x1/cp;
      phydro->w(IPR,j,i) = p0*pow(temp/Ts, cp/Rd);
      if (r <= 2.*a) {
        temp *= 1./(1. - qball + qball*pmicro->GetMassRatio(1));
        phydro->w(1,j,i) = qball;
      } else
        phydro->w(1,j,i) = 0.;
      phydro->w(IDN,j,i) = phydro->w(IPR,j,i)/(Rd*temp);
      phydro->w(IVX,j,i) = 0.;
      phydro->w(IVY,j,i) = 0.;
    }

  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je, ks, ke);
}
