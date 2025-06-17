# 3D Small-SFR: Steady-State

In this first case, a steady-state fully coupled neutronics, thermalhydraulics and thermomechanics simulation of a small SFR core is performed. The initial condition corresponds to the hot-zero-power configuration. Feedback from fuel temperature, cladding temperature, coolant density and core deformations are included. Moreover, the neutronics mesh is deformed to reflect the change of buckling due to the thermal expansions.

---

# Tasks

1. **Complete the regionsDict dictionary.** In regionDict, change the regionSolvers and mappings entries to correctly set-up the physics of the problem. In order to select the correct coupling fields, please refer to the field names table in https://foam-for-nuclear.gitlab.io/GeN-Foam/user_guide/coupling/
2. **Set up boundary conditions.** In particular:
   - Modify the 0/fluidRegion/U file to set the velocity BC and IC. Inlet velocity: 4.559114016 m/s. Use "slip" BC on the walls.
   - Modify the 0/fluidRegion/T file to set the temperature and IC. Inlet temperature: 668 K. Use "zeroGradient" elsewhere.
   - Modify the 0/neutroRegion/defaultFlux file to set the flux BC. For simplicity, use black boundary condition. 
3. **Set up the porous-media properties.** Modify the constant/fluidRegion/phaseProperties dictionary to assign the porous media properties to the fluid region. Follow the commented steps in the file to correctly:
   - For each cell-zone specify volume fraction and hydraulic diameter.
   - Define passive properties to model the thermal inertia of the structures.
   - Create a powerModel for the inner and outer core zones using the nuclearFuelPin model.
   - Assign drag and heat transfer coefficient models to all zones for both laminar and turbulent scenarios.
4. **Select neutronics calculation mode.** Modify the constant/neutroRegion/neutronicsProperties dictionary to select an eigenvalue-mode calculation.
5. **Specify simulation data.** Complete the system/controlDict file to perform a simulation of 200s, without adaptive time-stepping, with a delta-t = 1s, writing results every 100 s.
6. **Run case.** Run the simulation with GeN-Foam | tee log.GeN-Foam. To run in parallel do:
   - `decomposePar -allRegions`
   - `mpirun -np 4 GeN-Foam -parallel | tee log.GeN-Foam`