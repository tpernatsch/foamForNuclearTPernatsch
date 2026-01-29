## OFFBEAT Generic TRISO Fuel Particle Tutorial

This folder contains a tutorial for simulating a generic TRISO fuel particle using the OFFBEAT code. The tutorial demonstrates how OFFBEAT can be applied to study thermal and mechanical behavior in a 1D model of a TRISO fuel particle.

### TRISO Fuel Particle Model Overview

The generic TRISO fuel particle model consists of a 1D representation with double wedge geometry for spherical symmetry. It includes different layers representing the TRISO kernel, Buffer, IPyC (inner pyrolytic carbon), SiC (silicon carbide), and OPyC (outer pyrolytic carbon). Material properties and models are defined in the "materials" subdictionary inside the `solverDict`.

The geometry has been created using the sphereMaker.py script which can be adapted to other cases.

The tutorial employs the `smallStrainIncrementalUpdated` mechanics solver, which solves for the incremental displacement `DD`. This solver updates the mesh at every time step based on the converged solution from the previous time step.

For thermal analysis, a specific `trisoGap` boundary condition is used for heat transfer between the Buffer and IPyC layers. The tutorial also employs the TRISO version of the gapGas model to take into account the inner volume pressure during irradiation.

The volumetric heat generation rate (or vhgr) is set with the "timeDependentVhgr" heat source class. This allows to provide a list of time-dependent vhgr values directly in the `solverDict`.

### Instructions

The `Allrun` bash script can be used to create the mesh and automatically run the simulation. The `Allclean` script takes care of deleting all the unnecessary files from the case folder. 

The folder contains a python script called `create_mp4.py` that is able to convert a series of frames into an`.mp4` animation. In order to use the script, open ParaFoam and use the SaveAnimation option to save a certain number of frames (png files) into a folder called 'Frames'. Execute the Python script `create_mp4.py`. The script will convert the frames from the 'Frames' folder into an MP4 animation. Once the script finishes running, you will find the `animation.mp4` video in the 'Output' folder