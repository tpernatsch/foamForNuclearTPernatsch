# Description

In Case 3, we increase the complexity of our simulation and move to thermo-mechanics by activating the mechanical solver. The code now solves also for the radial displacement, following the 1.5D approximation in small strain, and the gap size changes over time due to thermal expansion. For this reason, we also include a gap model and related BC to include the effect of the changing gap size on the temperature distribution. Then, in Case 3a, we apply a single but important modification: starting from Case 3, we apply the isotropic cracking model so that the presence of cracks in the fuel pellet which affects the fuel stress state can be taken into account.

---

# Tasks for Case 3

1. **Create the folder.** Either copy `Case3` in your working directory or duplicate the `Case2` folder and rename it to `Case3` for this new simulation.
2. **Update Case3.** In particular:
   - Modify T Boundary Conditions or BC. Update the boundary conditions for `fuelOuter` and `cladInner` by selecting a `fuelRodGap` BC. This step is crucial as it adapts the thermal analysis to account for the changing gap size.
   - Activate mechanical analysis. Enable mechanical analysis in your simulation by choosing the `smallStrain` model as `mechanicsSolver`. Check the `mechanicsSolver` subdictionary options. Activate the `modifiedPlaneStrain` option in the `rheology` subdictionary in the `constant/solverDict` file so that the mechanical analysis follows the 1.5D approach. Use a `byMaterial` type for the `rheology` subdictionary.
   - Complete the mechanics BC in the `0/D` file. In particular, allow the `systemPressure` of 10MPa to act on the outer cladding surface.
   - Include `gapGas` model. Add a gap gas model by selecting the `Frapcon` model as `gapGas` in the `constant/solverDict` file.
   - Complete the `gapGas` file in the `0/` folder. Modify the file so that the gap contains only He and the initial pressure is 20 bars.
   - Add a `patchProbe` on the outer fuel patch. Complete the `fuelOuter` function in the `system/controlDict` file to register the `gapwidth` on the `fuelOuter` surface during the simulation.
3. **Run OFFBEAT Simulation.** Execute the OFFBEAT simulation* with your updated 'case3'. You can plot the results of the simulation using the provided `plot.py` Python script.
4. **Extract data and transform them for MATLAB.** Find and extract** the data stored at `postProcessing/fuelCenterline/0/T` (temperature at the center of the fuel over time), at `postProcessing/fuelOuter/0/gapWidth` (gapWidth at the outer fuel surface over time), at `postProcessing/fuelRadialProfile/##/line_T.xy` (temperature radial profile at ## time), and at `postProcessing/fuelRadialProfile/##/line_sigma.xy` (stress radial profile at ## time).
5. **Copy data in MATLAB grader.**

# Tasks for Case 3a

1. **Create the folder.** Duplicate the `Case3` folder and rename it to `Case3a` for this new simulation.
2. **Update Case3a.** In particular:
   - Activate isotropic cracking model. Enable the treatment of cracking in the fuel pellet by activating the `isotropicCracking` damage model (`damage` sub-dictionary within the material dict) in the fuel subdictionary within the materials section. Also, set a `nCracksMax` of 12, so that the model will be limited at a maximum crack number of 12 (indeed this is the standard version of the specific isotropic model used in OFFBEAT).
3. **Run OFFBEAT Simulation.** Execute the OFFBEAT simulation* with your updated 'case3a'. You can plot the results of the simulation using the provided `plot.py` Python script.
4. **Extract data and transform them for MATLAB.** Find and extract** the data stored at `postProcessing/fuelRadialProfile/##/line_sigma.xy` (stress radial profile at ## time).
5. **Copy only the radial and hoop stress for 168.04 days data in MATLAB grader.**

*NOTE: after you input the data and before running OFFBEAT, you need to set up the OpenFOAM geometry and create the mesh. For a 1D or 2D rod, the `blockMeshDict` file that is read by the `blockMesh` utility in order to create the mesh can be created with the `rodMaker.py` Python script provided within the case folder. For simplicity, all the required steps (from the creation of the geometry with `rodMaker` to running OFFBEAT) are included

 in the `Allrun` bash script that you can find in the folder and that can be conveniently run from terminal. If you need to start over for a new simulation, you can reset the folder using the `Allclean` bash script.

**NOTE: you can either extract the data manually following the instructions in class or use the `extractData` bash script for simplicity.