# 3D Small-SFR: Steady-State

In this case, a steady-state fully coupled neutronics,nthermalhydraulics and thermomechanics simulation of a small SFR core is performed. The initial condition corresponds to the hot-zero-power configuration. Feedback from fuel temperature, cladding temperature, coolant density and core deformations are included. The neutronics mesh is also deformed to reflect the effects of the change of geometrical buckling. 

---

# Tasks

1. **Complete the regionsDict dictionary.** In regionDict, change the regionSolvers and mappings entries to correctly set-up the physics of the problem. In order to select the correct coupling fields, please refer to the field names table in https://foam-for-nuclear.gitlab.io/GeN-Foam/user_guide/coupling/
2. **Run case.** Run the simulation with `GeN-Foam | tee log.GeN-Foam`. To run in parallel do:
   - `decomposePar -allRegions`
   - `mpirun -np 4 GeN-Foam -parallel | tee log.GeN-Foam`
3. **Plot Results** Run the "plotResults.py" file to study the convergence of k-eff and fuel and coolant temperatures. Compare the results to those in which core expansions were not included
4. **(Optional) TIghlty couple the physics** In system/regionsDict modify the *regionSolvers* entry as follows
```
regionSolvers
{
   Level_0
   {
      Level_1     picardLoop;
   }
   Level_1
   {
      subSolvers
      {
         fluidRegion             onePhase;
         neutroRegion            diffusionNeutronics;
         thermoMechanicalRegion  legacyThermoMechanics;
      }

      minResidual    1e-5;
      maxIterations  20;
   }
}
```
This will tightly coupled the physics together, performing Picard iterations until the convergence criterion is met. The compareResults.py can be used to compare also the tightly coupled to the loosely coupled results.