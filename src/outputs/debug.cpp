//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file debug.cpp
//  \brief writes debug output data, max and min quantities and their locations that are output
//         frequently in time to trace extreme values.

// C/C++ headers
#include <sstream>
#include <iostream>
#include <string>
#include <stdexcept>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <cfloat>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../globals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../hydro/hydro.hpp"
#include "../field/field.hpp"
#include "../mesh/mesh.hpp"
#include "outputs.hpp"

//----------------------------------------------------------------------------------------
// DebugOutput constructor
// destructor - not needed for this derived class

DebugOutput::DebugOutput(OutputParameters oparams)
  : OutputType(oparams)
{
}

//----------------------------------------------------------------------------------------
//! \fn void DebugOutput::WriteOutputFile()
//  \brief Writes a history file

void DebugOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin, bool flag)
{
  MeshBlock *pmb=pm->pblock;

  Real x1min[NHYDRO];
  Real x2min[NHYDRO];
  Real x3min[NHYDRO];
  Real x1max[NHYDRO];
  Real x2max[NHYDRO];
  Real x3max[NHYDRO];

  struct {
    Real value;
    Real coord;
  } lmin[NHYDRO], gmin[NHYDRO], lmax[NHYDRO], gmax[NHYDRO];

  for (int n = 0; n < NHYDRO; ++n) {
    lmin[n].value = FLT_MAX;
    lmax[n].value = FLT_MIN;
  }

  // Loop over MeshBlocks
  while (pmb != NULL) {
    Hydro *phydro = pmb->phydro;
    Coordinates *pcoord = pmb->pcoord;

    out_is=pmb->is; out_ie=pmb->ie;
    out_js=pmb->js; out_je=pmb->je;
    out_ks=pmb->ks; out_ke=pmb->ke;

    // ghost cells are included in the calculations
    out_is -= NGHOST; out_ie += NGHOST;
    if (out_js != out_je) {out_js -= NGHOST; out_je += NGHOST;}
    if (out_ks != out_ke) {out_ks -= NGHOST; out_ke += NGHOST;}

    // calculate maximum and minimum values over cells
    for (int n = 0; n < NHYDRO; ++n)
      for (int k=pmb->ks; k<=pmb->ke; ++k)
        for (int j=pmb->js; j<=pmb->je; ++j)
          for (int i=pmb->is; i<=pmb->ie; ++i) {
            if (phydro->w(n,k,j,i) < lmin[n].value) {
              lmin[n].value = phydro->w(n,k,j,i);
              x1min[n] = pcoord->x1v(i);
              x2min[n] = pcoord->x2v(j);
              x3min[n] = pcoord->x3v(k);
            }
            if (phydro->w(n,k,j,i) > lmax[n].value) {
              lmax[n].value = phydro->w(n,k,j,i);
              x1max[n] = pcoord->x1v(i);
              x2max[n] = pcoord->x2v(j);
              x3max[n] = pcoord->x3v(k);
            }
          }

    pmb = pmb->next;
  }

#ifdef MPI_PARALLEL
  // Because MPI_MINLOC/MPI_MAXLOC only allows two arguments.
  // I will call MPI three times, each gets one coordinate
  // x1 coordinate
  for (int n = 0; n < NHYDRO; ++n) lmin[n].coord = x1min[n];
  MPI_Reduce(lmin, gmin, NHYDRO, MPI_ATHENA_2REAL, MPI_MINLOC, 0, MPI_COMM_WORLD);
  for (int n = 0; n < NHYDRO; ++n) x1min[n] = gmin[n].coord;

  // x2 coordinate
  for (int n = 0; n < NHYDRO; ++n) lmin[n].coord = x2min[n];
  MPI_Reduce(lmin, gmin, NHYDRO, MPI_ATHENA_2REAL, MPI_MINLOC, 0, MPI_COMM_WORLD);
  for (int n = 0; n < NHYDRO; ++n) x2min[n] = gmin[n].coord;

  // x3 coordinate
  for (int n = 0; n < NHYDRO; ++n) lmin[n].coord = x3min[n];
  MPI_Reduce(lmin, gmin, NHYDRO, MPI_ATHENA_2REAL, MPI_MINLOC, 0, MPI_COMM_WORLD);
  for (int n = 0; n < NHYDRO; ++n) x3min[n] = gmin[n].coord;

  // now the root rank has all the extreme values and their coordinates
#endif

  // only the master rank writes the file
  // create filename: "file_basename" + ".hst".  There is no file number.
  if (Globals::my_rank == 0) {
    std::string fname;
    fname.assign(output_params.file_basename);
    fname.append(".debug");

    // open file for output
    FILE *pfile;
    std::stringstream msg;
    if((pfile = fopen(fname.c_str(),"a")) == NULL){
      msg << "### FATAL ERROR in function [OutputType::HistoryFile]" << std::endl
          << "Output file '" << fname << "' could not be opened";
      throw std::runtime_error(msg.str().c_str());
    }

    // If this is the first output, write header
    int iout = 1;
    if (output_params.file_number == 0)
      fprintf(pfile,"# Athena++ debug data\n"); // descriptor is first line

    fprintf(pfile,"# [%d]=time    ", iout++);
    fprintf(pfile,"[%d]=id        ", iout++);
    fprintf(pfile,"[%d]=dmin      ", iout++);
    fprintf(pfile,"[%d]=dmax      ", iout++);
    fprintf(pfile,"[%d]=x1min     ", iout++);
    fprintf(pfile,"[%d]=x2min     ", iout++);
    fprintf(pfile,"[%d]=x3min     ", iout++);
    fprintf(pfile,"[%d]=x1max     ", iout++);
    fprintf(pfile,"[%d]=x2max     ", iout++);
    fprintf(pfile,"[%d]=x3max     ", iout++);

    fprintf(pfile,"\n"); // terminate line

    // write debug variables
    for (int n = 0; n < NHYDRO; ++n) {
      if (n == 0)
        fprintf(pfile, output_params.data_format.c_str(), pm->time);
      else
        fprintf(pfile,"              ");
      fprintf(pfile, "%14d", n);
      fprintf(pfile, output_params.data_format.c_str(), gmin[n].value);
      fprintf(pfile, output_params.data_format.c_str(), gmax[n].value);
      fprintf(pfile, "%14d", x1min[n]);
      fprintf(pfile, "%14d", x2min[n]);
      fprintf(pfile, "%14d", x3min[n]);
      fprintf(pfile, "%14d", x1max[n]);
      fprintf(pfile, "%14d", x2max[n]);
      fprintf(pfile, "%14d", x3max[n]);
      fprintf(pfile,"\n");
    }
    fclose(pfile);
  }

  // increment counters, clean up
  output_params.file_number++;
  output_params.next_time += output_params.dt;
  pin->SetInteger(output_params.block_name, "file_number", output_params.file_number);
  pin->SetReal(output_params.block_name, "next_time", output_params.next_time);

  return;
}
