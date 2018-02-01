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
#include "../utils/utils.hpp"
#include "../math_funcs.hpp"

// EquationOfState constructor

EquationOfState::EquationOfState(MeshBlock *pmb, ParameterInput *pin)
{
  pmy_block_ = pmb;
  gamma_ = pin->GetReal("hydro","gamma");
  density_floor_  = pin->GetOrAddReal("hydro","dfloor",(1024*(FLT_MIN)));
  pressure_floor_ = pin->GetOrAddReal("hydro","pfloor",(1024*(FLT_MIN)));

  std::vector<std::string> str;

  // read eps, size should be NVAPOR
  if (NVAPOR > 0)
    SplitString(pin->GetOrAddString("hydro", "eps", ""), str);
  for (int n = 0; n < 1+2*NVAPOR; ++n) eps_[n] = 1.;
  for (int n = 0; n < std::min(NVAPOR, (int)str.size()); ++n) {
    eps_[1+n] = atof(str[n].c_str());     // gas
    eps_[1+2*n] = atof(str[n].c_str());   // condensate
  }

  // read rcv, size should be NVAPOR
  if (NVAPOR > 0)
    SplitString(pin->GetOrAddString("hydro", "rcv", ""), str);
  for (int n = 0; n < 1+2*NVAPOR; ++n) rcv_[n] = 1.;
  for (int n = 0; n < std::min(NVAPOR, (int)str.size()); ++n) {
    rcv_[1+n] = atof(str[n].c_str());     // gas
    rcv_[1+2*n] = atof(str[n].c_str());   // condensate
  }

  // read latent, size should be NVAPOR, unit is [J/kg]
  if (NVAPOR > 0)
    SplitString(pin->GetOrAddString("hydro", "latent", ""), str);
  for (int n = 0; n < 1+2*NVAPOR; ++n) latent_[n] = 0.;
  for (int n = 0; n < std::min(NVAPOR, (int)str.size()); ++n)
    latent_[1+2*n] = atof(str[n].c_str());
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

      Real& w_vx = prim(IVX,k,j,i);
      Real& w_vy = prim(IVY,k,j,i);
      Real& w_vz = prim(IVZ,k,j,i);
      Real& w_p  = prim(IEN,k,j,i);


      Real density = 0., sum = 0.;;
      for (int n = 0; n < ITR; ++n) {
        // apply density floor, without changing momentum or energy
        cons(n,k,j,i) = (cons(n,k,j,i) > density_floor_) ?  cons(n,k,j,i) : density_floor_;
        // record total density
        density += cons(n,k,j,i);
        prim(n,k,j,i) = cons(n,k,j,i)/eps_[n];
        sum += prim(n,k,j,i);
      }

      // density
      prim(IDN,k,j,i) = density;

      // volume mixing ratio
      for (int n = 1; n < ITR; ++n)
        prim(n,k,j,i) /= sum;

      // velocity
      Real di = 1.0/density;
      w_vx = u_m1*di;
      w_vy = u_m2*di;
      w_vz = u_m3*di;

      // internal energy
      Real KE = 0.5*di*(_sqr(u_m1) + _sqr(u_m2) + _sqr(u_m3));
      Real LE = 0.;
      for (int n = ICD; n < ICD + NVAPOR; ++n)
        LE += latent_[n]*cons(n,k,j,i);
      Real gmix = 1.;
      for (int n = 1; n < ITR; ++n)
        gmix += prim(n,k,j,i)*(rcv_[n] - 1.);
      w_p = gm1*(u_e - KE - LE)/gmix;

      // apply pressure floor, correct total energy
      u_e = (w_p > pressure_floor_) ?  u_e : ((pressure_floor_/gm1*gmix) + KE + LE);
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
      const Real& w_vx = prim(IVX,k,j,i);
      const Real& w_vy = prim(IVY,k,j,i);
      const Real& w_vz = prim(IVZ,k,j,i);
      const Real& w_p  = prim(IEN,k,j,i);

      // density
      Real qd = 1., sum = 0.;
      for (int n = 1; n < ITR; ++n) {
        cons(n,k,j,i) = prim(n,k,j,i)*eps_[n];
        qd -= prim(n,k,j,i);
        sum += cons(n,k,j,i);
      }
      cons(0,k,j,i) = qd;
      sum += qd;
      for (int n = 0; n < ITR; ++n)
        cons(n,k,j,i) = w_d*cons(n,k,j,i)/sum;

      // momentum
      u_m1 = w_vx*w_d;
      u_m2 = w_vy*w_d;
      u_m3 = w_vz*w_d;

      // total energy
      Real KE = 0.5*w_d*(_sqr(w_vx) + _sqr(w_vy) + _sqr(w_vz));
      Real LE = 0.;
      for (int n = ICD; n < ICD + NVAPOR; ++n)
        LE += latent_[n]*cons(n,k,j,i);
      Real gmix = 1.;
      for (int n = 1; n < ITR; ++n)
        gmix += prim(n,k,j,i)*(rcv_[n] - 1.);
      u_e = w_p/gm1*gmix + KE + LE;

      // set tracer density
      for (int n = ITR; n < ITR + NTRACER; ++n)
        cons(n,k,j,i) = w_d*prim(n,k,j,i);
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
  Real gmix = 1.;
  for (int n = 1; n < ITR; ++n)
    gmix += prim[n]*(rcv_[n] - 1.);
  return sqrt((gamma_ - 1. + gmix)/gmix*prim[IEN]/prim[IDN]);
}
