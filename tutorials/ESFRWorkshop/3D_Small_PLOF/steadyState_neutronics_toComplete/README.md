# 3D Small-SFR: Steady-State Stand-Alone Neutronics

In this first case, a steady-state stand-alone neutronics calculation is performed. Fuel and cladding temperature and coolant density are initialized to reflect hot-zero-power conditions

---

# Tasks

1. **Complete the regionsDict dictionary.** In regionDict, change the regionSolvers entry to correctly set-up the physics of the problem. 
2. **Set up boundary conditions.**  Modify the 0/neutroRegion/defaultFlux file to set the flux BC. For simplicity, use black boundary condition. 
3. **Select neutronics calculation mode.** Modify the constant/neutroRegion/neutronicsProperties dictionary to select an eigenvalue-mode calculation.
4. **Set up the reactorState dictionary.** Modify the 0/uniform file to set the initial guess for the eigenvalue and the reactor's power
5. **Select solution control specs.** In system/neutroRegion/fvSolution, change the neutronTransport entry to select 50 neutronics correctors with a minimum residual of 1e-6
6. **Run case.** Run the simulation with `GeN-Foam | tee log.GeN-Foam`
7. **Run again using SP3 instead** Change the consant/regionsDict file and select SP3Neutronics instead of diffusionNeutronics
8. **(Optional) Modify flux BC** Change the BC of defaultFlux and defaultFlux2 to the commented albedoSP3 BC in the file. By modifiying the "gamma" entry, one can change the albedo coefficient and see the change in k-eff 