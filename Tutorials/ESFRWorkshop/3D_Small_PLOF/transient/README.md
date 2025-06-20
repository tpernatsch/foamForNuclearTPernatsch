# 3D Small-SFR: PLOF Transient

Starting from the converged steady state, perform a transient simulation of a Proected Loss of Flow. This is achieved by decreasing the inlet velocity by 50% every 10 s until 25% of the nominal flow is reached. After 10 s from BOC, the reactor is scrammed by rapidly fully inserting all the CRs. 

---

# Tasks

1. **Copy the 200s folder from the steady state into the tranient folder.** The time-folder 200 previously obtained will be used as starting time folder for the transient simulation. After copying the file, please make sure to delete the 200/neutroRegion/polyMesh folder. This, in fact, corresponds to the deformed mesh file. If this mesh is used at restart:
   - The weights of the mesh-to-mesh mappings will be wrongly computed 
   - The displacement will be applied to an already-deformed mesh, duplicating the effect
2. **Modify the 200/fluidRegion/U boundary file.** The IC will be the non-uniform result from the steady-state. The inlet BC however needs to be modified to recreate a time-varying BC. This is achieved using the OpenFOAM *uniformFixedValue* BC as follows:
```
type              uniformFixedValue;
uniformValue      table
(
   (100     (0 0 4.559114))
   (150     (0 0 4.559114))
   (160     (0 0 2.279557))
   (170     (0 0 1.1397785))
   (2000    (0 0 1.1397785)) 
);
```
3. **Change the neutronics simulation mode.** Modify the constant/neutroRegion/neutronicsProperties dictionary to select a transient calculations
4. **Define the CRMove dictionary.** Complete the constant/neutroRegion/CRMove dictionary to simulate a full control rods insertion 10 s from BOC which is completed in 0.25 s. Assume the CR speed is 9 m/s.
5. **Specify simulation data.** Complete the system/controlDict file to perform a transient simulation of 300 s which starts from the last available time step, with adaptive time-stepping (limiting the Co number to 10 and the relative power variation to 0.025), with an initial delta-t = 1s, writing results every 10 s.
6. **Run case.** Run the simulation with GeN-Foam | tee log.GeN-Foam. To run in parallel do:
   - `decomposePar -allRegions`
   - `mpirun -np 4 GeN-Foam -parallel | tee log.GeN-Foam`
7. **Plot results** Plot the power and temperature evolution using the plotResults.py script
