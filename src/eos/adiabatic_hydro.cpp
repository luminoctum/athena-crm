//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file adiabatic_hydro.cpp
//  \brief implements functions in class EquationOfState for adiabatic hydrodynamics`

// C/C++ headers
#include <cmath>   // sqrt()
#include <cfloat>  // FLT_MIN

// Athena++ headers
#include "eos.hpp"
#include "../hydro/hydro.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../field/field.hpp"
#include "../math_funcs.hpp"
#include "../microphysics/microphysics.hpp"

// EquationOfState constructor

EquationOfState::EquationOfState(MeshBlock *pmb, ParameterInput *pin)
{
  pmy_block_ = pmb;
  gamma_ = pin->GetReal("hydro","gamma");
  density_floor_  = pin->GetOrAddReal("hydro","dfloor",(1024*(FLT_MIN)));
  pressure_floor_ = pin->GetOrAddReal("hydro","pfloor",(1024*(FLT_MIN)));
}

// destructor

EquationOfState::~EquationOfState()
{
}

//----------------------------------------------------------------------------------------
// \!fn void EquationOfState::ConservedToPrimitive(AthenaArray<Real> &cons,
//           const AthenaArray<Real> &prim_old, const FaceField &b,
//           AthenaArray<Real> &prim, AthenaArray<Real> &bcc, Coordinates *pco,
//           int is, int ie, int js, int je, int ks, int ke)
// \brief Converts conserved into primitive variables in adiabatic hydro.

void EquationOfState::ConservedToPrimitive(AthenaArray<Real> &cons,
  const AthenaArray<Real> &prim_old, const FaceField &b, AthenaArray<Real> &prim,
  AthenaArray<Real> &bcc, Coordinates *pco, int is, int ie, int js, int je, int ks, int ke)
{
  Real gm1 = GetGamma() - 1.0;
  Microphysics *pmicro = pmy_block_->pmicro;

  int nthreads = pmy_block_->pmy_mesh->GetNumMeshThreads();
#pragma omp parallel default(shared) num_threads(nthreads)
{
  for (int k=ks; k<=ke; ++k){
#pragma omp for schedule(dynamic)
  for (int j=js; j<=je; ++j){
#pragma simd
    for (int i=is; i<=ie; ++i){
      Real& u_m1 = cons(IM1,k,j,i);
      Real& u_m2 = cons(IM2,k,j,i);
      Real& u_m3 = cons(IM3,k,j,i);
      Real& u_e  = cons(IEN,k,j,i);

      Real& w_v1 = prim(IV1,k,j,i);
      Real& w_v2 = prim(IV2,k,j,i);
      Real& w_v3 = prim(IV3,k,j,i);
      Real& w_p  = prim(IPR,k,j,i);

      Real density = 0.;
      for (int n = 0; n < ITR; ++n) {
        // apply density floor, without changing momentum or energy
        cons(n,k,j,i) = (cons(n,k,j,i) > density_floor_) ?  cons(n,k,j,i) : density_floor_;
        // record total density
        density += cons(n,k,j,i);
      }

      // correct precipitation
      if (PRECIPITATION_ENABLED) {
        for (int n = ITR; n < ITR + NVAPOR; ++n)
          cons(n,k,j,i) = std::max(0., cons(n,k,j,i));
      }

      // total density
      prim(IDN,k,j,i) = density;
      Real di = 1./density;

      // mass mixing ratio
      for (int n = 1; n < ITR; ++n)
        prim(n,k,j,i) = cons(n,k,j,i)*di;

      // velocity
      w_v1 = u_m1*di;
      w_v2 = u_m2*di;
      w_v3 = u_m3*di;

      // internal energy
      Real KE = 0.5*di*(_sqr(u_m1) + _sqr(u_m2) + _sqr(u_m3));
      Real LE = 0., fsig = 1., feps = 1.;
      for (int n = ICD; n < ICD + NVAPOR; ++n) {
        LE += -pmicro->GetLatent(n)*cons(n,k,j,i);
        fsig += prim(n,k,j,i)*(pmicro->GetCvRatio(n) - 1.);
        feps -= prim(n,k,j,i);
      }
      for (int n = 1; n < ICD; ++n) {
        fsig += prim(n,k,j,i)*(pmicro->GetCvRatio(n) - 1.);
        feps += prim(n,k,j,i)*(1./pmicro->GetMassRatio(n) - 1.);
      }
      w_p = gm1*(u_e - KE - LE)*feps/fsig;

      // apply pressure floor, correct total energy
      u_e = (w_p > pressure_floor_) ?  u_e : (pressure_floor_/gm1*fsig/feps) + KE + LE;
      w_p = (w_p > pressure_floor_) ?  w_p : pressure_floor_;

      // set tracer mass mixing ratio
      for (int n = ITR; n < ITR + NTRACER; ++n)
        prim(n,k,j,i) = cons(n,k,j,i)*di;
    }
  }}
}

  return;
}


//----------------------------------------------------------------------------------------
// \!fn void EquationOfState::PrimitiveToConserved(const AthenaArray<Real> &prim,
//           const AthenaArray<Real> &bc, AthenaArray<Real> &cons, Coordinates *pco,
//           int is, int ie, int js, int je, int ks, int ke);
// \brief Converts primitive variables into conservative variables

void EquationOfState::PrimitiveToConserved(const AthenaArray<Real> &prim,
     const AthenaArray<Real> &bc, AthenaArray<Real> &cons, Coordinates *pco,
     int is, int ie, int js, int je, int ks, int ke)
{
  Real gm1 = GetGamma() - 1.0;
  Microphysics *pmicro = pmy_block_->pmicro;

  int nthreads = pmy_block_->pmy_mesh->GetNumMeshThreads();
#pragma omp parallel default(shared) num_threads(nthreads)
{
  for (int k=ks; k<=ke; ++k){
#pragma omp for schedule(dynamic)
  for (int j=js; j<=je; ++j){
#pragma simd
    for (int i=is; i<=ie; ++i){
      Real& u_m1 = cons(IM1,k,j,i);
      Real& u_m2 = cons(IM2,k,j,i);
      Real& u_m3 = cons(IM3,k,j,i);
      Real& u_e  = cons(IEN,k,j,i);

      const Real& w_d  = prim(IDN,k,j,i);
      const Real& w_v1 = prim(IV1,k,j,i);
      const Real& w_v2 = prim(IV2,k,j,i);
      const Real& w_v3 = prim(IV3,k,j,i);
      const Real& w_p  = prim(IPR,k,j,i);

      // density
      cons(IDN,k,j,i) = w_d;
      for (int n = 1; n < ITR; ++n) {
        cons(n,k,j,i) = prim(n,k,j,i)*w_d;
        cons(IDN,k,j,i) -= cons(n,k,j,i);
      }

      // momentum
      u_m1 = w_v1*w_d;
      u_m2 = w_v2*w_d;
      u_m3 = w_v3*w_d;

      // total energy
      Real KE = 0.5*w_d*(_sqr(w_v1) + _sqr(w_v2) + _sqr(w_v3));
      Real LE = 0., fsig = 1., feps = 1.;
      for (int n = ICD; n < ICD + NVAPOR; ++n) {
        LE += -pmicro->GetLatent(n)*cons(n,k,j,i);
        fsig += prim(n,k,j,i)*(pmicro->GetCvRatio(n) - 1.);
        feps -= prim(n,k,j,i);
      }
      for (int n = 1; n < ICD; ++n) {
        fsig += prim(n,k,j,i)*(pmicro->GetCvRatio(n) - 1.);
        feps += prim(n,k,j,i)*(1./pmicro->GetMassRatio(n) - 1.);
      }
      u_e = w_p/gm1*fsig/feps + KE + LE;

      // set tracer density
      for (int n = ITR; n < ITR + NTRACER; ++n)
        cons(n,k,j,i) = prim(n,k,j,i)*w_d;
    }
  }}
}
  return;
}

//----------------------------------------------------------------------------------------
// \!fn Real EquationOfState::SoundSpeed(Real prim[NHYDRO])
// \brief returns adiabatic sound speed given vector of primitive variables

Real EquationOfState::SoundSpeed(const Real prim[NHYDRO])
{
  Microphysics *pmicro = pmy_block_->pmicro;

  Real fsig = 1., feps = 1.;
  for (int n = ICD; n < ICD + NVAPOR; ++n) {
    fsig += prim[n]*(pmicro->GetCvRatio(n) - 1.);
    feps -= prim[n];
  }
  for (int n = 1; n < ICD; ++n) {
    fsig += prim[n]*(pmicro->GetCvRatio(n) - 1.);
    feps += prim[n]*(1./pmicro->GetMassRatio(n) - 1.);
  }

  return sqrt((1. + (gamma_ - 1)*feps/fsig)*prim[IPR]/prim[IDN]);
}
