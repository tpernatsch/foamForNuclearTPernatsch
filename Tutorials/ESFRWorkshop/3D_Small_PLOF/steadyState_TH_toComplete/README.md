# 3D Small-SFR: Steady-State

In this case, a steady-state stand-alone thermalhydraulics simulation of a small SFR core is performed. The power density obtained from the neutronics simulation is used to calculate the coolant temperature distribution 

---

# Tasks

1. **Complete the regionsDict dictionary.** In regionDict, change the regionSolvers entry to correctly set-up the physics of the problem. 
2. **Set up boundary conditions.** In particular:
   - Modify the 0/fluidRegion/U file to set the velocity BC and IC. Inlet velocity: 4.559114016 m/s. Use "slip" BC on the walls.
   - Modify the 0/fluidRegion/T file to set the temperature and IC. Inlet temperature: 668 K. Use "zeroGradient" elsewhere.
3. **Set up the porous-media properties.** Modify the constant/fluidRegion/phaseProperties dictionary to assign the porous media properties to the fluid region. Follow the commented steps in the file to correctly:
   - For each cell-zone specify volume fraction and hydraulic diameter.
   - Define passive properties to model the thermal inertia of the structures.
   - Create a powerModel for the inner and outer core zones using the nuclearFuelPin model.
   - Assign drag and heat transfer coefficient models to all zones for both laminar and turbulent scenarios.
4. **Specify PIMPLE loop specs.** In the system/fluidRegion/fvSolution file, modify the PIMPLE entry to perform PIMPLE loops with two correctors and 6 outer correctors
6. **Run case.** Run the simulation with GeN-Foam | tee log.GeN-Foam. To run in parallel do:
   - `decomposePar -region fluidRegion`
   - `mpirun -np 4 GeN-Foam -parallel | tee log.GeN-Foam`