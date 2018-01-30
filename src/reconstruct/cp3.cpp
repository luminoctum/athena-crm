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
inline Real interp_cp3(Real phim1, Real phi, Real phip1) {
  return 1./3.*phim1 + 5./6.*phi - 1./6.*phip1;
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
      ql(n,i) = interp_cp3(q(n,k,j,i),q(n,k,j,i-1),q(n,k,j,i-2));
      qr(n,i) = interp_cp3(q(n,k,j,i-1),q(n,k,j,i),q(n,k,j,i+1));
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
      ql(n,i) = interp_cp3(q(n,k,j,i),q(n,k,j-1,i),q(n,k,j-2,i));
      qr(n,i) = interp_cp3(q(n,k,j-1,i),q(n,k,j,i),q(n,k,j+1,i));
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
      ql(n,i) = interp_cp3(q(n,k,j,i),q(n,k-1,j,i),q(n,k-2,j,i));
      qr(n,i) = interp_cp3(q(n,k-1,j,i),q(n,k,j,i),q(n,k+1,j,i));
    }
  }

  return;
}
