## BiMaterial Hot Thin Plate Tutorial

This folder contains a tutorial for simulating a biMaterial hot thin plate using the OFFBEAT code. 
The tutorial aims to guide users in understanding how to model the deformation of such body.

### BiMaterial Hot Thin Plate Overview

The biMaterial hot thin plate consists of two materials with different properties, listed in the `materials` subdictionary inside `solverDict`. 
When exposed to temperatures above the reference temperature of 300K, the plate experiences bending due to the differential thermal expansion of the materials.

### Case Setup

To model the biMaterial hot thin plate efficiently, we use symmetry boundary conditions to simulate only a quarter of the plate. This approach reduces computational effort while providing meaningful insights into the bending behavior.

For the analysis, we activate only the mechanics and thermal solvers. We use the `smallStrain` mechanics solver, which solves the momentum equation for the total displacement field D always in the original mesh (i.e. the mesh is not updated).
Both sides of the plate are modeled as linear-elastic materials.

### Instructions

The `Allrun` bash script can be used to create the mesh and automatically run the simulation. The `Allclean` script takes care of deleting all the unnecessary files from the case folder. 

The folder contains a python script called `create_mp4.py` that is able to convert a series of frames into an`.mp4` animation. In order to use the script, open ParaFoam and use the SaveAnimation option to save a certain number of frames (png files) into a folder called 'Frames'. Execute the Python script `create_mp4.py`. The script will convert the frames from the 'Frames' folder into an MP4 animation. Once the script finishes running, you will find the `bend.mp4` video in the 'Output' folder