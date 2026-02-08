# Description

This exercise involves three sub-cases, each a variation of the original Case1. These variations are designed to explore different aspects of thermal behavior in a nuclear fuel rod simulation without the complexity of mechanics and nuclear-specific phenomena.

## Sub-case Description:
- **Case1a**: Similar to Case1, but with a different constant value for the fuel thermal conductivity of $k_f = 2 W/(m\cdot K)$
- **Case1b**: Similar to Case1, but with a different gap conductance, equal to 2000 $W/(m^2 \cdot K)$
- **Case1c**: Similar to Case1, but the fuel in this scenario has a central hole with a radius of 0.75mm, altering its thermal profile.

---

# Tasks for each sub-case

1. **Create the folder.** Copy the completed `Case1`, change name (e.g. `Case1a`) and clean the files using the `Allclean` bash script to prepare for the new simulation.
2. **Update the Case folder.** Use the details from the provided above*.
3. **Run OFFBEAT Simulation.** Execute the OFFBEAT simulation* with your updated `case1`. You can plot the results of the simulation using the provided `plot.py` Python script.
4. **Extract data and transform them for MATLAB.** Find and extract** the data stored at postProcessing/fuelCenterline/0/T (temperature at the center of the fuel over time) and at `postProcessing/fuelRadialProfile/##/line_T.xy` (temperature radial profile at ## time).
5. **Copy data in MATLAB grader.**

*NOTE: `Case1c` has a central hole. For this reason, the probes must be slightly changed compared to the reference case. Use the code snippet below in the functions section of the `system/controlDict` file.
**NOTE: after you input the data and before running OFFBEAT, you need to set up the OpenFOAM geometry and create the mesh. For a 1D or 2D rod, the `blockMeshDict` file that is read by the `blockMesh` utility in order to create the mesh can be created with the `rodMaker.py` Python script provided within the case folder. For simplicity, all the required steps (from the creation of the geometry with `rodMaker` to running OFFBEAT) are included in the `Allrun` bash script that you can find in the folder and that can be conveniently run from terminal. If you need to start over for a new simulation, you can reset the folder using the `Allclean` bash script.
*** NOTE: you can either extract the data manually following the instructions in class or use the `extractData` bash script for simplicity.

## Code snippet 
```text
functions 
{ 
    fuelCenterline
    {
        type            patchProbes;
        libs            ("libsampling.so");
        
        writeControl    timeStep;
        writeInterval   1;
        
        probeLocations
        (
            (0.00075 0.0 0.055)
        );
        patchName fuelInner;
        
        fields  (T);
    }

    fuelRadialProfile
    {
        type            sets;
        libs            ("libsampling.so");
        
        writeControl    writeTime;
        interpolationScheme cell;
        setFormat       raw;
        
        sets
        (
            line
            {
                type            lineUniform;
                axis            x;
                start           (0.00075 0.0 0.055);
                end             (0.0045 0.0 0.055);
                nPoints         30;
            }
        );  

        fields          (T);      
    }
}
