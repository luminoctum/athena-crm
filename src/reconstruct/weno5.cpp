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
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"

// WENO 5 interpolation
inline Real interp_weno5(Real phim2, Real phim1, Real phi, Real phip1, Real phip2) {
  Real p0 = (1.0/3.0)*phim2 - (7.0/6.0)*phim1 + (11.0/6.0)*phi;
  Real p1 = (-1.0/6.0) * phim1 + (5.0/6.0)*phi + (1.0/3.0)*phip1;
  Real p2 = (1.0/3.0) * phi + (5.0/6.0)*phip1 - (1.0/6.0)*phip2;

  Real beta2 = (13.0/12.0 * (phi - 2.0 * phip1 + phip2)*(phi - 2.0 * phip1 + phip2)
                        + 0.25 * (3.0 * phi - 4.0 * phip1 + phip2)*(3.0 * phi - 4.0 * phip1 + phip2));
  Real beta1 = (13.0/12.0 * (phim1 - 2.0 * phi + phip1)*(phim1 - 2.0 * phi + phip1)
                        + 0.25 * (phim1 - phip1)*(phim1 - phip1));
  Real beta0 = (13.0/12.0 * (phim2 - 2.0 * phim1 + phi)*(phim2 - 2.0 * phim1 + phi)
                        + 0.25 * (phim2 - 4.0 * phim1 + 3.0 * phi)*(phim2 - 4.0 * phim1 + 3.0 * phi));

  Real alpha0 = 0.1/((beta0 + 1e-10) * (beta0 + 1e-10));
  Real alpha1 = 0.6/((beta1 + 1e-10) * (beta1 + 1e-10));
  Real alpha2 = 0.3/((beta2 + 1e-10) * (beta2 + 1e-10));

  Real alpha_sum_inv = 1.0/(alpha0 + alpha1 + alpha2);

  Real w0 = alpha0 * alpha_sum_inv;
  Real w1 = alpha1 * alpha_sum_inv;
  Real w2 = alpha2 * alpha_sum_inv;

  return w0 * p0 + w1 * p1 + w2 * p2;
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
