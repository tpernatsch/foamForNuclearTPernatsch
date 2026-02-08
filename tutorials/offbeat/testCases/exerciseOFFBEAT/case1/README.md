# Description

In the next few exercises, we're going to look at a short rodlet made of 10 pellets that undergoes a simplified one-year irradiation. We'll start with the basics, neglecting nuclear-specific phenomena such as relocation, densification or creep-down and considering only thermal analysis: think of the fuel inside the rodlet as just a hot cylinder that gives off heat to the surrounding Zircaloy cladding through a gap with a fixed width and heat conductance.

---

# Data

## Geometry:
- The pellet has no central hole
- Fuel outer radius 4.5 $mm$
- Gap size 60 $\mu m$
- Radial clad width 700 $\mu m$
- Fuel column made of 10 pellets, each  tall 11 $mm$
- The top plenum is  tall 10 $mm$
- We neglect bottom and top cap

## Discretization:
- 30 cells in the fuel radially, 1 axially
- 10 cells in the clad radially, 1 axially

## Temperature data:
- Initial temperature $T_0=293K$
- Cladding outer temperature $T_{c0}=573K$
- Fixed gap conductance $\alpha_{gap}$ or $h_{gap}=5000 W/(m^2 \cdot K)$  

## Linear heat rate (lhgr) during 1 year of operation:
- Linear ramp-up from 0 to 40 $kW/m$ in 0.04 days (approx. 1 hr)
- Power-hold at 40 $kW/m$ until 364 days
- Linear ramp-down from 40 to 0 $kW/m$ in 0.04 days (approx. 1 hr)
- Cool-down phase until 365 days

## Material properties:
- Fuel thermal conductivity $k_f = 3 W/(m \cdot K)$
- Clad thermal conductivity $k_c = 21 W/(m \cdot K)$

---

# Tasks

1. **Copy Case1 in your working directory.**
2. **Update Case1.** Use the details from the provided data sheet to update your simulation setup. In particular:
   - Complete dimensions and mesh settings in the `rodDict` file.
   - Complete initial and boundary conditions in `0/T` file.
   - Complete material properties and `lhgr` in constant/`solverDict` file
3. **Run OFFBEAT Simulation.** Execute the OFFBEAT simulation* with your updated '`case1`'. You can plot the results of the simulation using the provided `plot.py` Python script.
4. **Extract data and transform them for MATLAB.** Find and extract** the data stored at `postProcessing/fuelCenterline/0/T` (temperature at the center of the fuel over time) and at `postProcessing/fuelRadialProfile/##/line_T.xy` (temperature radial profile at ## time).
5. **Copy data in MATLAB grader.**

*NOTE: after you input the data and before running OFFBEAT, you need to set up the OpenFOAM geometry and create the mesh. For a 1D or 2D rod, the `blockMeshDict` file that is read by the blockMesh utility in order to create the mesh can be created with the `rodMaker.py` Python script provided within the case folder. For simpicity, all the required steps (from the creation of the geometry with `rodMaker` to running OFFBEAT) are included in the `Allrun` bash script that you can find in the folder and that can be conveniently run from terminal. If you need to start over for a new simulation, you can reset the folder using the `Allclean` bash script.

** NOTE: you can either extract the data manually following the instructionts in class or use the `extractData` bash script for simplicity.
