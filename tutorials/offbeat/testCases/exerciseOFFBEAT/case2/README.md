# Description

In Case 2, we start incorporating nuclear-specific phenomena into our thermal analysis. While maintaining the same settings as in Case1, a key difference here is that thermal conductivity is no longer constant but varies as a function of temperature. Additionally, as for the UO2 fuel conductivity is also significantly influenced by burnup, we activate the burnup module which integrates the volumetric power density (times the density) over time.

---

# Tasks

1. **Create the folder.** Either copy `Case2` in your working directory or duplicate the `Case1` folder and rename it to `Case2` for this new simulation.
2. **Update Case2.** In particular:
   - Update conductivity model by removing in the `constant/solverDict` file the `conductivity` sub-dicts. This action will prompt the code to use the default conductivity correlations for the two materials.
   - Introduce a burnup solver in the `solverDict` file (`burnup` sub-dict), using the `fromPower` model `type`. This model calculates burnup based on the integrated volumetric power density over time.
   - Modify the probes by adding `Bu` (burnup) to the `fuelCenterline` probe in `system/controlDict`.
3. **Run OFFBEAT Simulation.** Execute the OFFBEAT simulation* with your updated 'case1'. You can plot the results of the simulation using the provided `plot.py` Python script.
4. **Extract data and transform them for MATLAB.** Find and extract** the data stored at `postProcessing/fuelCenterline/0/T` (temperature at the center of the fuel over time) and at `postProcessing/fuelRadialProfile/##/line_T.xy` (temperature radial profile at ## time).
5. **Copy data in MATLAB grader.**

*NOTE: after you input the data and before running OFFBEAT, you need to set up the OpenFOAM geometry and create the mesh. For a 1D or 2D rod, the `blockMeshDict` file that is read by the `blockMesh` utility in order to create the mesh can be created with the `rodMaker.py` Python script provided within the case folder. For simplicity, all the required steps (from the creation of the geometry with `rodMaker` to running OFFBEAT) are included in the `Allrun` bash script that you can find in the folder and that can be conveniently run from terminal. If you need to start over for a new simulation, you can reset the folder using the `Allclean` bash script.

**NOTE: you can either extract the data manually following the instructions in class or use the `extractData` bash script for simplicity.