## BiMaterial Hot Cylinder Tutorial

This folder contains a tutorial for simulating a biMaterial hot cylinder using the OFFBEAT code. 
The tutorial aims to guide users in understanding how to model the deformation of such body.

### BiMaterial Hot Cylinder Overview

The biMaterial hot cylinder consists of two materials with different properties, listed in the `materials` subdictionary inside `solverDict`. 
The cylinder is exposed to a temperature and pressure gradient from the inner (800 K, 1e7 Pa) to the outer boundary (300 K, 1e5 Pa).

### Case Setup

To model the biMaterial hot cylinder efficiently, we use cylic boundary conditions to simulate only a portion of the cylinder. This approach reduces computational effort while providing meaningful insights into the thermo-mechanical behavior.

For the analysis, we activate only the mechanics and thermal solvers. We use the `smallStrain` mechanics solver, which solves the momentum equation for the total displacement field D always in the original mesh (i.e. the mesh is not updated).
Both sides of the cylinder are modeled as linear-elastic materials.

### Instructions

The `Allrun` bash script can be used to create the mesh and automatically run the simulation. The `Allclean` script takes care of deleting all the unnecessary files from the case folder. 