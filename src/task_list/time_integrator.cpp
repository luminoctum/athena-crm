//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file time_integrator.cpp
//  \brief derived class for time integrator task list.  Can create task lists for one
//  of many different time integrators (e.g. van Leer, RK2, RK3, etc.)

// C/C++ headers
#include <iostream>   // endl
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ classes headers
#include "task_list.hpp"
#include "../athena.hpp"
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../field/field.hpp"
#include "../bvals/bvals.hpp"
#include "../eos/eos.hpp"
#include "../hydro/srcterms/hydro_srcterms.hpp"
#include "../microphysics/microphysics.hpp"

//----------------------------------------------------------------------------------------
//  TimeIntegratorTaskList constructor

TimeIntegratorTaskList::TimeIntegratorTaskList(ParameterInput *pin, Mesh *pm)
  : TaskList(pm)
{
  // First, set weights for each step of time-integration algorithm.  Each step is
  //    U^{2} = a*U^0 + b*U^1 + c*dt*Div(F), where U^0 and U^1 are previous steps 
  // a,b=(1-a),and c are weights that are different for each step and each integrator
  // These are stored as: time_int_wght1 = a, time_int_wght2 = b, time_int_wght3 = c

  integrator = pin->GetOrAddString("time","integrator","rk3");

  // second-order van Leer integrator (Gardiner & Stone, NewA 14, 139 2009)
  if (integrator == "vl2") {
    nsub_steps = 2;
    step_wghts[0].a = 1.0;
    step_wghts[0].b = 0.0;
    step_wghts[0].c = 0.5;

    step_wghts[1].a = 1.0;
    step_wghts[1].b = 0.0;
    step_wghts[1].c = 1.0;
  } else if (integrator == "rk2") {
    nsub_steps = 2;
    step_wghts[0].a = 1.0;
    step_wghts[0].b = 0.0;
    step_wghts[0].c = 1.0;

    step_wghts[1].a = 0.5;
    step_wghts[1].b = 0.5;
    step_wghts[1].c = 0.5;
  } else if (integrator == "rk3") {
    nsub_steps = 3;
    step_wghts[0].a = 1.0;
    step_wghts[0].b = 0.0;
    step_wghts[0].c = 1.0;

    step_wghts[1].a = 0.75;
    step_wghts[1].b = 0.25;
    step_wghts[1].c = 0.25;

    step_wghts[2].a = 1./3.;
    step_wghts[2].b = 2./3.;
    step_wghts[2].c = 2./3.;
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in CreateTimeIntegrator" << std::endl
        << "integrator=" << integrator << " not valid time integrator" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  // Now assemble list of tasks for each step of time integrator
  {using namespace HydroIntegratorTaskNames;
    AddTimeIntegratorTask(START_ALLRECV,NONE);

    // compute hydro fluxes, integrate hydro variables
    AddTimeIntegratorTask(CALC_HYDFLX,START_ALLRECV);
    if(pm->multilevel==true) { // SMR or AMR
      AddTimeIntegratorTask(SEND_HYDFLX,CALC_HYDFLX);
      AddTimeIntegratorTask(RECV_HYDFLX,CALC_HYDFLX);
      AddTimeIntegratorTask(INT_HYD, RECV_HYDFLX);
    } else {
      AddTimeIntegratorTask(INT_HYD, CALC_HYDFLX);
    }
    AddTimeIntegratorTask(SRCTERM_HYD,INT_HYD);
    AddTimeIntegratorTask(SEND_HYD,SRCTERM_HYD);
    AddTimeIntegratorTask(RECV_HYD,START_ALLRECV);

    // compute MHD fluxes, integrate field
    if (MAGNETIC_FIELDS_ENABLED) { // MHD
      AddTimeIntegratorTask(CALC_FLDFLX,CALC_HYDFLX);
      AddTimeIntegratorTask(SEND_FLDFLX,CALC_FLDFLX);
      AddTimeIntegratorTask(RECV_FLDFLX,SEND_FLDFLX);
      AddTimeIntegratorTask(INT_FLD, RECV_FLDFLX);
      AddTimeIntegratorTask(SEND_FLD,INT_FLD);
      AddTimeIntegratorTask(RECV_FLD,START_ALLRECV);
    }

    // prolongate, compute new primitives
    if (MAGNETIC_FIELDS_ENABLED) { // MHD
      if(pm->multilevel==true) { // SMR or AMR
        AddTimeIntegratorTask(PROLONG, (SEND_HYD|RECV_HYD|SEND_FLD|RECV_FLD));
        AddTimeIntegratorTask(CON2PRIM,PROLONG);
      } else {
        AddTimeIntegratorTask(CON2PRIM,(INT_HYD|RECV_HYD|INT_FLD|RECV_FLD));
      }
    } else {  // HYDRO
      if(pm->multilevel==true) { // SMR or AMR
        AddTimeIntegratorTask(PROLONG,(SEND_HYD|RECV_HYD));
        AddTimeIntegratorTask(CON2PRIM,PROLONG);
      } else {
        AddTimeIntegratorTask(CON2PRIM,(INT_HYD|RECV_HYD));
      }
    }

    // microphysics
    AddTimeIntegratorTask(MICROPHY,CON2PRIM);

    // everything else
    //AddTimeIntegratorTask(PHY_BVAL,CON2PRIM);
    AddTimeIntegratorTask(PHY_BVAL,MICROPHY);
    AddTimeIntegratorTask(USERWORK,PHY_BVAL);
    AddTimeIntegratorTask(NEW_DT,USERWORK);
    if(pm->adaptive==true) {
      AddTimeIntegratorTask(AMR_FLAG,USERWORK);
      AddTimeIntegratorTask(CLEAR_ALLRECV,AMR_FLAG);
    } else {
      AddTimeIntegratorTask(CLEAR_ALLRECV,NEW_DT);
    }

  } // end of using namespace block
}

//----------------------------------------------------------------------------------------//! \fn
//  \brief Sets id and dependency for "ntask" member of task_list_ array, then iterates
//  value of ntask.  

void TimeIntegratorTaskList::AddTimeIntegratorTask(uint64_t id, uint64_t dep)
{
  task_list_[ntasks].task_id=id;
  task_list_[ntasks].dependency=dep;

  using namespace HydroIntegratorTaskNames;
  switch((id)) {
    case (START_ALLRECV):
      task_list_[ntasks].TaskFunc= 
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::StartAllReceive);
      break;
    case (CLEAR_ALLRECV):
      task_list_[ntasks].TaskFunc= 
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::ClearAllReceive);
      break;

    case (CALC_HYDFLX):
      task_list_[ntasks].TaskFunc= 
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::CalculateFluxes);
      break;
    case (CALC_FLDFLX):
      task_list_[ntasks].TaskFunc= 
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::CalculateEMF);
      break;

    case (SEND_HYDFLX):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::FluxCorrectSend);
      break;
    case (SEND_FLDFLX):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::EMFCorrectSend);
      break;

    case (RECV_HYDFLX):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::FluxCorrectReceive);
      break;
    case (RECV_FLDFLX):
      task_list_[ntasks].TaskFunc= 
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::EMFCorrectReceive);
      break;

    case (INT_HYD):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::HydroIntegrate);
      break;
    case (INT_FLD):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::FieldIntegrate);
      break;

    case (SRCTERM_HYD):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::HydroSourceTerms);
      break;

    case (SEND_HYD):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::HydroSend);
      break;
    case (SEND_FLD):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::FieldSend);
      break;

    case (RECV_HYD):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::HydroReceive);
      break;
    case (RECV_FLD):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::FieldReceive);
      break;

    case (PROLONG):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::Prolongation);
      break;
    case (CON2PRIM):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::Primitives);
      break;
    case (MICROPHY):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::ApplyMicrophysics);
      break;
    case (PHY_BVAL):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::PhysicalBoundary);
      break;
    case (USERWORK):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::UserWork);
      break;
    case (NEW_DT):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::NewBlockTimeStep);
      break;
    case (AMR_FLAG):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::CheckRefinement);
      break;

    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in AddTask" << std::endl
          << "Invalid Task "<< id << " is specified" << std::endl;
      throw std::runtime_error(msg.str().c_str());
  }
  ntasks++;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn
//  \brief

//----------------------------------------------------------------------------------------
// Functions to start/end MPI communication

enum TaskStatus TimeIntegratorTaskList::StartAllReceive(MeshBlock *pmb, int step)
{
  pmb->pbval->StartReceivingAll();
  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskList::ClearAllReceive(MeshBlock *pmb, int step)
{
  pmb->pbval->ClearBoundaryAll();
  return TASK_SUCCESS;
}

//----------------------------------------------------------------------------------------
// Functions to calculates fluxes

enum TaskStatus TimeIntegratorTaskList::CalculateFluxes(MeshBlock *pmb, int step)
{
  Hydro *phydro=pmb->phydro;
  Field *pfield=pmb->pfield;

  if (step == 1) {
    // use first order reconstruction if vl2 integrator and first step
    if (integrator == "vl2")
      phydro->CalculateFluxes(phydro->w,  pfield->b,  pfield->bcc, true);
    else
      phydro->CalculateFluxes(phydro->w,  pfield->b,  pfield->bcc, false);
    return TASK_NEXT;
  } else {
    phydro->CalculateFluxes(phydro->w1, pfield->b1, pfield->bcc1, false);
    return TASK_NEXT;
  }

  return TASK_FAIL;
}

enum TaskStatus TimeIntegratorTaskList::CalculateEMF(MeshBlock *pmb, int step)
{
  if(step == 1) {
    pmb->pfield->ComputeCornerE(pmb->phydro->w,  pmb->pfield->bcc);
    return TASK_NEXT;
  }

  if(step == 2) {
    pmb->pfield->ComputeCornerE(pmb->phydro->w1, pmb->pfield->bcc1);
    return TASK_NEXT;
  } 

  return TASK_FAIL;
}

//----------------------------------------------------------------------------------------
// Functions to communicate fluxes between MeshBlocks for flux correction step with AMR

enum TaskStatus TimeIntegratorTaskList::FluxCorrectSend(MeshBlock *pmb, int step)
{
  pmb->pbval->SendFluxCorrection();
  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskList::EMFCorrectSend(MeshBlock *pmb, int step)
{
  pmb->pbval->SendEMFCorrection();
  return TASK_SUCCESS;
}

//----------------------------------------------------------------------------------------
// Functions to receive fluxes between MeshBlocks

enum TaskStatus TimeIntegratorTaskList::FluxCorrectReceive(MeshBlock *pmb, int step)
{
  if(pmb->pbval->ReceiveFluxCorrection() == true) {
    return TASK_NEXT;
  } else {
    return TASK_FAIL;
  }
}

enum TaskStatus TimeIntegratorTaskList::EMFCorrectReceive(MeshBlock *pmb, int step)
{
  if(pmb->pbval->ReceiveEMFCorrection() == true) {
    return TASK_NEXT;
  } else {
    return TASK_FAIL;
  }
}

//----------------------------------------------------------------------------------------
// Functions to integrate conserved variables

enum TaskStatus TimeIntegratorTaskList::HydroIntegrate(MeshBlock *pmb, int step)
{
  Hydro *ph=pmb->phydro;
  Field *pf=pmb->pfield;

  if(step == 1) {
    ph->AddFluxDivergenceToAverage(ph->u,ph->u,ph->w,pf->bcc,step_wghts[step-1],ph->u1);
    return TASK_NEXT;
  } else if (step == nsub_steps) {
    ph->AddFluxDivergenceToAverage(ph->u,ph->u1,ph->w1,pf->bcc1,step_wghts[step-1],ph->u);
    return TASK_NEXT;
  } else {
    ph->AddFluxDivergenceToAverage(ph->u,ph->u1,ph->w1,pf->bcc1,step_wghts[step-1],ph->u1);
    return TASK_NEXT;
  }

  return TASK_FAIL;
}

enum TaskStatus TimeIntegratorTaskList::FieldIntegrate(MeshBlock *pmb, int step)
{
  if(step == 1) {
    pmb->pfield->CT(pmb->pfield->b, pmb->pfield->b, step_wghts[0], pmb->pfield->b1);
    return TASK_NEXT;
  }

  if((step == 2) && (integrator == "vl2")) {
    pmb->pfield->CT(pmb->pfield->b, pmb->pfield->b, step_wghts[1], pmb->pfield->b);
    return TASK_NEXT;
  }

  if((step == 2) && (integrator == "rk2")) {
    pmb->pfield->CT(pmb->pfield->b, pmb->pfield->b1, step_wghts[1], pmb->pfield->b);
    return TASK_NEXT;
  }

  return TASK_FAIL;
}

//----------------------------------------------------------------------------------------
// Functions to add source terms

enum TaskStatus TimeIntegratorTaskList::HydroSourceTerms(MeshBlock *pmb, int step)
{
  Hydro *ph=pmb->phydro;
  Field *pf=pmb->pfield;

  // return if there are no source terms to be added
  if (ph->psrc->hydro_sourceterms_defined == false) return TASK_NEXT;

  Real dt = (step_wghts[(step-1)].c)*(pmb->pmy_mesh->dt);
  Real time;

  if (step == 1) {
    time=pmb->pmy_mesh->time;
    ph->psrc->AddHydroSourceTerms(time,dt,ph->flux,ph->w,pf->bcc,ph->u1);
  } else if (step == nsub_steps) {
    time = pmb->pmy_mesh->time + (step_wghts[step-2].c)*pmb->pmy_mesh->dt;
    ph->psrc->AddHydroSourceTerms(time,dt,ph->flux,ph->w1,pf->bcc1,ph->u);
  } else {
    time = pmb->pmy_mesh->time + (step_wghts[step-2].c)*pmb->pmy_mesh->dt;
    ph->psrc->AddHydroSourceTerms(time,dt,ph->flux,ph->w1,pf->bcc1,ph->u1);
  }

  return TASK_NEXT;
}

//----------------------------------------------------------------------------------------
// Functions to communicate conserved variables between MeshBlocks

enum TaskStatus TimeIntegratorTaskList::HydroSend(MeshBlock *pmb, int step)
{
  if(step != nsub_steps) {
    pmb->pbval->SendHydroBoundaryBuffers(pmb->phydro->u1, true);
  } else { // last step
    pmb->pbval->SendHydroBoundaryBuffers(pmb->phydro->u, true);
  }

  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskList::FieldSend(MeshBlock *pmb, int step)
{
  if(step == 1) {
    pmb->pbval->SendFieldBoundaryBuffers(pmb->pfield->b1);
  } else if(step == 2) {
    pmb->pbval->SendFieldBoundaryBuffers(pmb->pfield->b);
  } else {
    return TASK_FAIL;
  }
  return TASK_SUCCESS;
}

//----------------------------------------------------------------------------------------
// Functions to receive conserved variables between MeshBlocks

enum TaskStatus TimeIntegratorTaskList::HydroReceive(MeshBlock *pmb, int step)
{
  bool ret;
  if(step != nsub_steps) {
    ret=pmb->pbval->ReceiveHydroBoundaryBuffers(pmb->phydro->u1);
  } else { // last step
    ret=pmb->pbval->ReceiveHydroBoundaryBuffers(pmb->phydro->u);
  }
  if(ret==true) {
    return TASK_SUCCESS;
  } else {
    return TASK_FAIL;
  }
}

enum TaskStatus TimeIntegratorTaskList::FieldReceive(MeshBlock *pmb, int step)
{
  bool ret;
  if(step == 1) {
    ret=pmb->pbval->ReceiveFieldBoundaryBuffers(pmb->pfield->b1);
  } else if(step == 2) {
    ret=pmb->pbval->ReceiveFieldBoundaryBuffers(pmb->pfield->b);
  } else {
    return TASK_FAIL;
  }
  if(ret==true) {
    return TASK_SUCCESS;
  } else {
    return TASK_FAIL;
  }
}

//----------------------------------------------------------------------------------------
// Functions for everything else

enum TaskStatus TimeIntegratorTaskList::Prolongation(MeshBlock *pmb, int step)
{
  Hydro *phydro=pmb->phydro;
  Field *pfield=pmb->pfield;
  BoundaryValues *pbval=pmb->pbval;
  Real dt;

  if(step == 1) {
    dt = (step_wghts[(step-1)].c)*(pmb->pmy_mesh->dt);
    pbval->ProlongateBoundaries(phydro->w1, phydro->u1, pfield->b1, pfield->bcc1,
                                pmb->pmy_mesh->time+dt, dt);
  } else if(step == 2) {
    dt=pmb->pmy_mesh->dt;
    pbval->ProlongateBoundaries(phydro->w,  phydro->u,  pfield->b,  pfield->bcc,
                                pmb->pmy_mesh->time+dt, dt);
  } else {
    return TASK_FAIL;
  }
  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskList::Primitives(MeshBlock *pmb, int step)
{
  Hydro *phydro=pmb->phydro;
  Field *pfield=pmb->pfield;
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  if(pmb->nblevel[1][1][0]!=-1) is-=NGHOST;
  if(pmb->nblevel[1][1][2]!=-1) ie+=NGHOST;
  if(pmb->nblevel[1][0][1]!=-1) js-=NGHOST;
  if(pmb->nblevel[1][2][1]!=-1) je+=NGHOST;
  if(pmb->nblevel[0][1][1]!=-1) ks-=NGHOST;
  if(pmb->nblevel[2][1][1]!=-1) ke+=NGHOST;

  if(step != nsub_steps) {
    pmb->peos->ConservedToPrimitive(phydro->u1, phydro->w, pfield->b1,
                                    phydro->w1, pfield->bcc1, pmb->pcoord,
                                    is, ie, js, je, ks, ke);
  } else {  // last step
    pmb->peos->ConservedToPrimitive(phydro->u, phydro->w1, pfield->b,
                                    phydro->w, pfield->bcc, pmb->pcoord,
                                    is, ie, js, je, ks, ke);
  }

  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskList::ApplyMicrophysics(MeshBlock *pmb, int step)
{
  Microphysics *pmicro = pmb->pmicro;
  Hydro *ph = pmb->phydro;
  Field *pf = pmb->pfield;
  // return if there are no source terms to be added
  if (pmicro->ncycle == 0) return TASK_NEXT;

  // only do microphysics at last rk step
  if (step < nsub_steps) return TASK_NEXT;

  // do microphysics every xx step
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  if (pmb->pmy_mesh->ncycle % pmicro->ncycle == 0) {
    pmicro->SaturationAdjustment(ph->w,is,ie,js,je,ks,ke);
    if (PRECIPITATION_ENABLED) {
      Real dt = pmb->pmy_mesh->dt*pmicro->ncycle;
      pmicro->Precipitation(ph->w,ph->u,dt,is,ie,js,je,ks,ke);
      pmicro->Evaporation(ph->w,ph->u,dt);
    }
    pmb->peos->PrimitiveToConserved(ph->w, pf->bcc, ph->u, pmb->pcoord,
      is, ie, js, je, ks, ke);
  }

  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskList::PhysicalBoundary(MeshBlock *pmb, int step)
{
  Hydro *phydro=pmb->phydro;
  Field *pfield=pmb->pfield;
  BoundaryValues *pbval=pmb->pbval;
  Real dt;
  if(step != nsub_steps) {
    dt = (step_wghts[(step-1)].c)*(pmb->pmy_mesh->dt);
    pbval->ApplyPhysicalBoundaries(phydro->w1, phydro->u1, pfield->b1, pfield->bcc1,
                                   pmb->pmy_mesh->time+dt, dt);
  } else {  // last step
    dt=pmb->pmy_mesh->dt;
    pbval->ApplyPhysicalBoundaries(phydro->w,  phydro->u,  pfield->b,  pfield->bcc,
                                   pmb->pmy_mesh->time+dt, dt);
  }

  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskList::UserWork(MeshBlock *pmb, int step)
{
  if (step != nsub_steps) return TASK_SUCCESS; // only do on last sub-step

  pmb->UserWorkInLoop();
  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskList::NewBlockTimeStep(MeshBlock *pmb, int step)
{
  if (step != nsub_steps) return TASK_SUCCESS; // only do on last sub-step

  pmb->phydro->NewBlockTimeStep();
  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskList::CheckRefinement(MeshBlock *pmb, int step)
{
  if (step != nsub_steps) return TASK_SUCCESS; // only do on last sub-step

  pmb->pmr->CheckRefinementCondition();
  return TASK_SUCCESS;
}
