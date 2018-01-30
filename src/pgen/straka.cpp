//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file straka.cpp
//  \brief Problem generator for RT instabilty.
//
// Note the gravitational acceleration is hardwired to be 0.1. Density difference is
// hardwired to be 2.0 in 2D, and is set by the input parameter <problem>/rhoh in 3D
// (default value is 3.0). This reproduces 2D results of Liska & Wendroff, 3D results of
// Dimonte et al.
// 
// FOR 2D HYDRO:
// Problem domain should be -1/6 < x < 1/6; -0.5 < y < 0.5 with gamma=1.4 to match Liska
// & Wendroff. Interface is at y=0; perturbation added to Vy. Gravity acts in y-dirn.
// Special reflecting boundary conditions added in x2 to improve hydrostatic eqm
// (prevents launching of weak waves) Atwood number A=(d2-d1)/(d2+d1)=1/3. Options:
//     iprob = 1  -- Perturb V2 using single mode
//     iprob != 1 -- Perturb V2 using multiple mode
//
// FOR 3D:
// Problem domain should be -.05 < x < .05; -.05 < y < .05, -.1 < z < .1, gamma=5/3 to
// match Dimonte et al.  Interface is at z=0; perturbation added to Vz. Gravity acts in
// z-dirn. Special reflecting boundary conditions added in x3.  A=1/2.  Options:
//     iprob = 1 -- Perturb V3 using single mode
//     iprob = 2 -- Perturb V3 using multiple mode
//     iprob = 3 -- B rotated by "angle" at interface, multimode perturbation
//
// REFERENCE: R. Liska & B. Wendroff, SIAM J. Sci. Comput., 25, 995 (2003)
//========================================================================================

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../utils/utils.hpp"
#include "../globals.hpp"
#include "../math_funcs.hpp"

Real grav, p0, Ts, Rd, cp, K;

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  grav = -phydro->psrc->GetG1();
  p0 = pin->GetReal("problem", "p0");
  Ts = pin->GetReal("problem", "Ts");
  Rd = pin->GetReal("problem", "Rd");
  cp = pin->GetReal("problem", "cp");
  K = pin->GetReal("problem", "K");

  Real xc = pin->GetReal("problem", "xc");
  Real xr = pin->GetReal("problem", "xr");
  Real zc = pin->GetReal("problem", "zc");
  Real zr = pin->GetReal("problem", "zr");
  Real dT = pin->GetReal("problem", "dT");

  for (int j = js; j <= je; ++j)
    for (int i = is; i <= ie; ++i) {
      Real x1 = pcoord->x1v(i);
      Real x2 = pcoord->x2v(j);
      Real L = sqrt(_sqr((x2 - xc)/xr) + _sqr((x1 - zc)/zr));
      Real temp = Ts - grav*x1/cp;
      phydro->w(IPR,j,i) = p0*pow(temp/Ts, cp/Rd);
      if (L <= 1.)
        temp += dT*(cos(M_PI*L) + 1.)/2.;
      phydro->w(IDN,j,i) = phydro->w(IPR,j,i)/(Rd*temp);
      phydro->w(IVX,j,i) = 0.;
      phydro->w(IVY,j,i) = 0.;
    }

  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je, ks, ke);
}
