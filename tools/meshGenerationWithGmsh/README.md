# Mesh generation with GMSH

This folder contains three examples of mesh generation for a SFR core using [gmsh](https://gmsh.info/).

To obtain each mesh one should:
1. enter the folder
2. type `gmsh esfr.main` in the terminal
3. wait (a long time)
4. mesh the geometry once the window is unfrozed
5. save the `.msh` file that can then be converted to OpenFOAM using `gmshToFoam`

Different core geometries can be obtained by changing parameters in `coreGeometry.geo`. 

These are just supposed to be examples, not general-purpose tools for geometry creation. Comments are limited.
