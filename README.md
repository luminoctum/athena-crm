# Athena/CRM: Athena++ cloud resolving model
## Minimum recipe to change the original Athena++ to a cloud resolving model
1. add rk3 integrator  
File changed:
* src/task_list/time_integrator.cpp

2. modify reconstruction and add weno3, weno5, cp3, cp5 methods  
File changed:
* src/reconstruct/plm.cpp
* src/reconstruct/reconstruction.hpp
* src/reconstruct/reconstruction.cpp
File added:
* src/reconstruct/cp3.pp
* src/reconstruct/cp5.pp
* src/reconstruct/weno3.pp
* src/reconstruct/weno5.pp

3. add netcdf output option  
File changed:
* src/outputs/outputs.hpp
* src/outputs/outputs.cpp
File added:
* src/outputs/netcdf.cpp

4. add two static pressure variables in hydro  
File changed:
* src/hydro/hydro.hpp
* src/hydro/hydro.cpp

5. add pressure decomposition subroutines  
File changed:
* src/hydro/calculate_fluxes.cpp
File added:
* src/hydro/decompose_assemble.cpp

6. add LMARS riemann solver  
File added:
* src/hydro/rsolvers/hydro/lmars.cpp

7. add utility files:  
File changed:
* src/utils/utils.hpp
File aded:  
* src/utils/utils.cpp

8. add math function files:  
File added:
* src/math_funcs.hpp

9. change definitions/outputs accordingly:  
File changed:
* configure.py
* Makefile.in
* src/defs.hpp.in
* src/utils/show_config.cpp

10. add straka and robert test case:  
File added:
* src/pgen/straka.cpp
* src/pgen/robert.cpp

11. add Coriolis acceleration: 681ef4d  
File changed:
* src/hydro/srcterms/hydro_srcterms.hpp
* src/hydro/srcterms/hydro_srcterms.cpp
File added:
* src/hydro/srcterms/coriolis_acc.cpp

12. add tracer advection (only works for adiabatic_hydro now)  
File changed:
* configure.py
* src/defs.hpp.in
* src/athena.hpp
* src/eos/adiabatic_hydro.cpp
* src/hydro/hydro.hpp
* src/hydro/hydro.cpp
* src/hydro/calculate_fluxes.cpp
* src/outputs/outputs.cpp

13. add vapor/cloud (only works for adiabatic_hydro now)   
File changed:
* configure.py
* src/defs.hpp.in
* src/athena.hpp
* src/eos/eos.hpp
* src/eos/adiabatic_hydro.cpp
* src/hydro/rsolvers/hydro/lmars.cpp
* src/outputs/outputs.cpp

14. add microphysics  
File changed:
* src/Makefile.in
* src/task_list/task_list.hpp
* src/task_list/time_integrator.cpp
* src/eos/eos.hpp (this file is reverted to the file in step 13.)
* src/eos/adiabatic_hydro.cpp
* src/mesh/mesh.hpp
* src/mesh/meshblock.cpp
* src/outputs/outputs.cpp
File added:
* src/microphysics/microphysics.hpp 
* src/microphysics/microphysics.cpp 
* src/microphysics/saturation_adjust_jupiter.cpp 
* src/pgen/jad.cpp
* src/utils/logging_info.hpp

## Known bugs (but didn't have time to correct it)
1. flip_across_pole_hydro
Only works with standard hydro
2. src/output/output.cpp:
I don't know whether I should use pdata = pnew->pnext or pdata = pdata->pnext
