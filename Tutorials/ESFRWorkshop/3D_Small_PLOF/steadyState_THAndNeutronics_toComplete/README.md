# 3D Small-SFR: Steady-State

In this case, a steady-state fully coupled neutronics and thermalhydraulics simulation of a small SFR core is performed. The initial condition corresponds to the hot-zero-power configuration. Feedback from fuel temperature, cladding temperature, coolant density are included. 

---

# Tasks

1. **Complete the regionsDict dictionary.** In regionDict, change the regionSolvers and mappings entries to correctly set-up the physics of the problem. In order to select the correct coupling fields, please refer to the field names table in https://foam-for-nuclear.gitlab.io/GeN-Foam/user_guide/coupling/
2. **Activate remove baffles options.** In system/controlDict, activate the removeBaffles option for the fluidRegion. This is necessary to correctly perform the mapping between fluid and neutroncis regions
3. **Run case.** Run the simulation with GeN-Foam | tee log.GeN-Foam. To run in parallel do:
   - `decomposePar -allRegions`
   - `mpirun -np 4 GeN-Foam -parallel | tee log.GeN-Foam`
4. **Plot Results** Run the "plotResults.py" file to study the convergence of k-eff and fuel and coolant temperatures