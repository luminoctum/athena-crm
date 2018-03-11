//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file plm.cpp
//  \brief  piecewise linear reconstruction

// Athena++ headers
#include "reconstruction.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../math_funcs.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
//#include "../coordinates/coordinates.hpp"

// WENO 5 interpolation
inline void interp_weno5(Real phim2, Real phim1, Real phi, Real phip1, Real phip2, Real *phiL, Real *phiR) {
  Real p0R = (1.0/3.0)*phim2 - (7.0/6.0)*phim1 + (11.0/6.0)*phi;
  Real p1R = (-1.0/6.0)*phim1 + (5.0/6.0)*phi + (1.0/3.0)*phip1;
  Real p2R = (1.0/3.0)*phi + (5.0/6.0)*phip1 - (1.0/6.0)*phip2;

  Real p0L = (1.0/3.0)*phip2 - (7.0/6.0)*phip1 + (11.0/6.0)*phi;
  Real p1L = (-1.0/6.0)*phip1 + (5.0/6.0)*phi + (1.0/3.0)*phim1;
  Real p2L = (1.0/3.0)*phi + (5.0/6.0)*phim1 - (1.0/6.0)*phim2;

  Real beta0 = 13./12.*_sqr(phim2 - 2.*phim1 + phi) + 0.25*_sqr(phim2 - 4.*phim1 + 3.*phi);
  Real beta1 = 13./12.*_sqr(phim1 - 2.*phi + phip1) + 0.25*_sqr(phim1 - phip1);
  Real beta2 = 13./12.*_sqr(phi - 2.*phip1 + phip2) + 0.25*_sqr(3.*phi - 4.*phip1 + phip2);

  Real alpha0 = .1/_sqr(beta0 + 1e-6);
  Real alpha1 = .6/_sqr(beta1 + 1e-6);
  Real alpha2 = .3/_sqr(beta2 + 1e-6);

  *phiR = (alpha0*p0R + alpha1*p1R + alpha2*p2R)/(alpha0 + alpha1 + alpha2);
  *phiL = (alpha2*p0L/3. + alpha1*p1L + alpha0*p2L*3.)/(alpha2/3. + alpha1 + alpha0*3.);
};

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::ReconstructionFuncX1()
//  \brief 

void Reconstruction::HighResFuncX1(const int k, const int j,
  const int il, const int iu,
  const AthenaArray<Real> &q, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &ql, AthenaArray<Real> &qr)
{
  for (int n=0; n<NHYDRO; ++n)
    for (int i=il-1; i<=iu; ++i) {
      interp_weno5(q(n,k,j,i-2),q(n,k,j,i-1),q(n,k,j,i),q(n,k,j,i+1),q(n,k,j,i+2),&qr(n,i),&ql(n,i+1));
    }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::ReconstructionFuncX2()
//  \brief 

void Reconstruction::HighResFuncX2(const int i, const int k,
  const int jl, const int ju,
  const AthenaArray<Real> &q, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &ql, AthenaArray<Real> &qr)
{
  for (int n=0; n<NHYDRO; ++n)
    for (int j=jl-1; j<=ju; ++j){
      interp_weno5(q(n,k,j-2,i),q(n,k,j-1,i),q(n,k,j,i),q(n,k,j+1,i),q(n,k,j+2,i),&qr(n,j),&ql(n,j+1));
    }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::ReconstructionFuncX3()
//  \brief 

void Reconstruction::HighResFuncX3(const int k, const int j,
  const int il, const int iu,
  const AthenaArray<Real> &q, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &ql, AthenaArray<Real> &qr)
{
  for (int n=0; n<NHYDRO; ++n) {
#pragma simd
    for (int i=il; i<=iu; ++i){
      //ql(n,i) = interp_weno5(q(n,k-3,j,i),q(n,k-2,j,i),q(n,k-1,j,i),q(n,k,j,i),q(n,k+1,j,i));
      //qr(n,i) = interp_weno5(q(n,k+2,j,i),q(n,k+1,j,i),q(n,k,j,i),q(n,k-1,j,i),q(n,k-2,j,i));
    }
  }

  return;
}
