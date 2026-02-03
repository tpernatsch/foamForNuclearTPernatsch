# Description

In Case 5, we increase further the complexity of our simulation by adding the phenomenon of cladding creep, moving beyond the purely elastic behavior observed in the earlier cases. Due to creep, the Zircaloy cladding accumulates plastic (or permanent) strain when exposed to moderately high temperatures and stresses over extended periods of time. The impact of neutron fast flux, which is known to enhance the creep phenomenon, is now a crucial factor to consider.

Due to differential pressure exerted on the cladding's inner and outer surfaces - typically higher on the outside, especially at lower fission gas release fractions - cladding creep results in the material "creeping" down towards the fuel pellets. This movement, known as cladding creep-down, is a key contributor to the closure of the gap between the fuel and the cladding, thus influencing the overall behavior of the fuel rod under operational conditions.

---

# Tasks

1. **Create the folder.** Either copy `Case5` in your working directory or duplicate the `Case4` folder and rename it to `Case5` for this new simulation.
2. **Update Case5.** In particular:
   - Change mechanical constitutive law for the cladding. Modify the `constitutiveLaw` to `misesPlasticCreep` in the cladding subdictionary within the materials section of the `constant/solverDict` file.
   - Copmlete the options in `constitutiveLaw` setting a high constant yield stress (e.g. 400 MPa). Within this subdictionary create an additional dict for `creep` and set the `ZircaloyLimback` as `type`, set a relax factor of 0.9 and define a constant yield stress of 400MPa.
   - Add a `fastFlux` subdictionary  in `solverDict` file, selecting `timeDependentAxialProfile` as `type`.
   - Complete entries in `fastFlux` subdictionary (see Code Snippets section). Assuming a proportional relationship between fast flux and power, linearly ramp up the fast flux from 0 to 1e13 n/cm²/s in 0.04 days (approximately 1 hour) and maintain it constant until day 364. Then, decrease the fast flux to 0 in 0.04 days and keep it at 0 until 365 days.
3. **Run OFFBEAT Simulation.** Execute the OFFBEAT simulation* with your updated 'case5'. You can plot the results of the simulation using the provided `plot.py` Python script.
4. **Extract data and transform them for MATLAB.** Find and extract** the data stored at `postProcessing/fuelCenterline/0/T` (temperature at the center of the fuel over time), at `postProcessing/fuelOuter/0/gapWidth` (gapWidth at the outer fuel surface over time).
5. **Copy data in MATLAB grader.**

*NOTE: after you input the data and before running OFFBEAT, you need to set up the OpenFOAM geometry and create the mesh. For a 1D or 2D rod, the `blockMeshDict` file that is read by the `blockMesh` utility in order to create the mesh can be created with the `rodMaker.py` Python script provided within the case folder. For simplicity, all the required steps (from the creation of the geometry with `rodMaker` to running OFFBEAT) are included in the `Allrun` bash script that you can find in the folder and that can be conveniently run from terminal. If you need to start over for a new simulation, you can reset the folder using the `Allclean` bash script.

**NOTE: you can either extract the data manually following the instructions in class or use the `extractData` bash script for simplicity.

## Code Snippet #1
```text
fastFlux
{
    type timeDependentAxialProfile;
    timePoints  ( 0  0.4    364    364.04 365 );
    fastFlux    ( 0  1e13  1e13    0      0);
    
    axialProfile
    {
        type flat;
    }

    materials ( fuel cladding );
}
```