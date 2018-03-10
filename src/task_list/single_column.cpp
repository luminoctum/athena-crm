// C/C++ headers
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

// Athena++ class headers
#include "task_list.hpp"
#include "../athena.hpp"
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../coordinates/coordinates.hpp"
#include "../microphysics/microphysics.hpp"

// global variables
using namespace SingleColumnTaskNames;

SingleColumnTaskList::SingleColumnTaskList(ParameterInput *pin, Mesh *pm):
  TaskList(pm)
{
  nsub_steps = 1;

  AddTask(CONST_AD, NONE);
}

void SingleColumnTaskList::AddTask(uint64_t id, uint64_t dep)
{
  task_list_[ntasks].task_id = id;
  task_list_[ntasks].dependency = dep;

  switch (id) {
    case (CONST_AD):
      task_list_[ntasks].TaskFunc = 
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&SingleColumnTaskList::ConstructAdiabat);
      break;

    default:
      std::stringstream msg; 
      msg << "### FATAL ERROR in AddTask" << std::endl
          << "Invalid Task " << id << " is specified" << std::endl;
      throw std::runtime_error(msg.str().c_str());
  }
  ntasks++;
  return;
}

enum TaskStatus SingleColumnTaskList::ConstructAdiabat(MeshBlock *pmb, int step)
{
  /*Microphysics *pmicro = pmb->pmicro;
  Hydro *phydro = pmb->phydro;
  Coordiante *pcoord = pmb->pcoord;

  int ks = pmb->ks, js = pmb->js, is = pmb->is;
  int ke = pmb->ke, je = pmb->je, ie = pmb->ie;
  int i0 = _locate(pcoord->x1v.data(), 0., pcoord->x1v.GetSize());

  pmicro->MoistAdiabat(phydro->w, T0, P0, grav, ks, js, i0, is, ie)*/

  return TASK_SUCCESS;
}
