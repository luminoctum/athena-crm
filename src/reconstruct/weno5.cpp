//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file weno5.cpp
//  \brief  WENO5 interpolation

// Athena++ headers
#include "reconstruction.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../math_funcs.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
//#include "../coordinates/coordinates.hpp"

// WENO 5 interpolation
inline Real interp_weno5(Real phim2, Real phim1, Real phi, Real phip1, Real phip2) {
  Real p0 = (1./3.)*phim2 - (7./6.)*phim1 + (11./6.)*phi;
  Real p1 = (-1./6.)*phim1 + (5./6.)*phi + (1./3.)*phip1;
  Real p2 = (1./3.)*phi + (5./6.)*phip1 - (1./6.)*phip2;

  Real beta0 = 13./12.*_sqr(phim2 - 2.*phim1 + phi) + .25*_sqr(phim2 - 4.*phim1 + 3.*phi);
  Real beta1 = 13./12.*_sqr(phim1 - 2.*phi + phip1) + .25*_sqr(phim1 - phip1);
  Real beta2 = 13./12.*_sqr(phi - 2.*phip1 + phip2) + .25*_sqr(3.*phi - 4.*phip1 + phip2);

  Real alpha0 = .1/_sqr(beta0 + 1e-6);
  Real alpha1 = .6/_sqr(beta1 + 1e-6);
  Real alpha2 = .3/_sqr(beta2 + 1e-6);

  return (alpha0*p0 + alpha1*p1 + alpha2*p2)/(alpha0 + alpha1 + alpha2);
};

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::ReconstructionFuncX1()
//  \brief 

void Reconstruction::HighResFuncX1(const int k, const int j,
  const int il, const int iu,
  const AthenaArray<Real> &q, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &ql, AthenaArray<Real> &qr)
{
  for (int n=0; n<NHYDRO; ++n) {
#pragma simd
    for (int i=il; i<=iu; ++i){
      ql(n,i) = interp_weno5(q(n,k,j,i-3),q(n,k,j,i-2),q(n,k,j,i-1),q(n,k,j,i),q(n,k,j,i+1));
      qr(n,i) = interp_weno5(q(n,k,j,i+2),q(n,k,j,i+1),q(n,k,j,i),q(n,k,j,i-1),q(n,k,j,i-2));
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::ReconstructionFuncX2()
//  \brief 

void Reconstruction::HighResFuncX2(const int k, const int j,
  const int il, const int iu,
  const AthenaArray<Real> &q, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &ql, AthenaArray<Real> &qr)
{
  for (int n=0; n<NHYDRO; ++n) {
#pragma simd
    for (int i=il; i<=iu; ++i){
      ql(n,i) = interp_weno5(q(n,k,j-3,i),q(n,k,j-2,i),q(n,k,j-1,i),q(n,k,j,i),q(n,k,j+1,i));
      qr(n,i) = interp_weno5(q(n,k,j+2,i),q(n,k,j+1,i),q(n,k,j,i),q(n,k,j-1,i),q(n,k,j-2,i));
    }
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
      ql(n,i) = interp_weno5(q(n,k-3,j,i),q(n,k-2,j,i),q(n,k-1,j,i),q(n,k,j,i),q(n,k+1,j,i));
      qr(n,i) = interp_weno5(q(n,k+2,j,i),q(n,k+1,j,i),q(n,k,j,i),q(n,k-1,j,i),q(n,k-2,j,i));
    }
  }

  return;
}
