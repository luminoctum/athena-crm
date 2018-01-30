//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file reconstruction.cpp
//  \brief 

// Athena++ headers
#include <iostream>
#include "reconstruction.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp" 
#include "../mesh/mesh.hpp"

// constructor

Reconstruction::Reconstruction(MeshBlock *pmb, ParameterInput *pin)
{
  pmy_block_ = pmb;
  int nxmax = std::max(pmb->block_size.nx1,
    std::max(pmb->block_size.nx2, pmb->block_size.nx3));
  nxmax += 2*NGHOST;
  c1_.NewAthenaArray(nxmax);
  c2_.NewAthenaArray(nxmax);
  c3_.NewAthenaArray(nxmax);
  c4_.NewAthenaArray(nxmax);
  c5_.NewAthenaArray(nxmax);
  c6_.NewAthenaArray(nxmax);
  dd.NewAthenaArray(nxmax);
  dph.NewAthenaArray(nxmax);
}

// destructor

Reconstruction::~Reconstruction()
{
  c1_.DeleteAthenaArray();
  c2_.DeleteAthenaArray();
  c3_.DeleteAthenaArray();
  c4_.DeleteAthenaArray();
  c5_.DeleteAthenaArray();
  c6_.DeleteAthenaArray();
  dd.DeleteAthenaArray();
  dph.DeleteAthenaArray();
}
