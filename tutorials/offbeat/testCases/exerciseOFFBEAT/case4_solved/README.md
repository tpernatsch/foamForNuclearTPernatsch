# Description

Case 4 is identical to Case 3, but it introduces two crucial nuclear fuel phenomena: densification and relocation. Densification refers to the closing of porosity in nuclear fuel leading to pellet shrinkage and gap opening. Relocation is the outward movement of fuel fragments post-cracking. This phenomenon is considered the major contributor to gap closing.

---

# Tasks

1. **Create the folder.** Either copy `Case4` in your working directory or duplicate the `Case3a` folder and rename it to `Case4` for this new simulation.
2. **Update Case4.** In particular:
   - Introduce densification model. Identify the `densification` sub-dictionary in the fuel subdictionary within the materials section of the `constant/solverDict` file. Select the `empirical` model as type.
   - Complete the densification model settings. Introduce a `densificationDensityChange` of 1% (given as a percent) and a `densificationTimeConstant` of 1000 MWd/tU.
   - Introduce relocation model. Identify the `relocation` subdictionary in the fuel subdictionary within the materials section of the `constant/solverDict` file. Select the `UO2Frapcon` model as `type`.
   - Complete the relocation model settings. Introduce the `GapCold` and `DiamCold` parameters (given in m) and the name of the `outerPatch` of the fuel.
3. **Run OFFBEAT Simulation.** Execute the OFFBEAT simulation* with your updated 'case4'. You can plot the results of the simulation using the provided `plot.py` Python script.
4. **Extract data and transform them for MATLAB.** Find and extract** the data stored at `postProcessing/fuelCenterline/0/T` (temperature at the center of the fuel over time), at `postProcessing/fuelOuter/0/gapWidth` (gapWidth at the outer fuel surface over time).
5. **Copy data in MATLAB grader.**

*NOTE: after you input the data and before running OFFBEAT, you need to set up the OpenFOAM geometry and create the mesh. For a 1D or 2D rod, the `blockMeshDict` file that is read by the `blockMesh` utility in order to create the mesh can be created with the `rodMaker.py` Python script provided within the case folder. For simplicity, all the required steps (from the creation of the geometry with `rodMaker` to running OFFBEAT) are included in the `Allrun` bash script that you can find in the folder and that can be conveniently run from terminal. If you need to start over for a new simulation, you can reset the folder using the `Allclean` bash script.

**NOTE: you can either extract the data manually following the instructions in class or use the `extractData` bash script for simplicity.