## OFFBEAT Pseudo-PUZRY Tutorial

This folder contains a pseudo-PUZRY tutorial for the OFFBEAT code. The tutorial is designed to help users understand how to use OFFBEAT for studying separate-effect tests such as burst tests.

### PUZRY Series Overview

The PUZRY series comprises short cladding segments that undergo different internal pressure increases at varying rates and are kept at different temperatures. These tests were conducted to study high-temperature creep and helped the derivation of creep correlation for burst analysis.

### OFFBEAT Case Description

The OFFBEAT case consists of a short wedge or axisymmetric model that represents half of the cladding segment (i.e. bottom part). Symmetry boundary conditions are applied to the top surface.

Only mechanics and thermal analysis are activated for this OFFBEAT simulation. The `smallStrainIncrementalUpdated` mechanics solver is used. This means that the momentum equation is solved for the incremental displacement field DD, for which internal field and boundary conditions should be set in the `0/` folder. Also the updated version of the solver means that the mesh is moved at the beginning of every time-step, using the converged solution obtained in the previous time-step.

The rheological behavior is modeled using the `misesPlasticCreep` model, with `ZircaloyLimbackLoca` as the creep model. Notably, the primary and irradiation creep have been turned off for consistency with simulations conducted by other authors with the MOOSE platform.

### Overstrain Failure Criterion

An overstrain failure criterion is incorporated, which checks the maximum strain in the outer and inner patches.

### Instructions

The `Allrun` bash script can be used to create the mesh and automatically run the simulation. The `Allclean` script takes care of deleting all the unnecessary files from the case folder. 

The folder contains a python script called `create_mp4.py` that is able to convert a series of frames into an`.mp4` animation. In order to use the script, open ParaFoam and use the SaveAnimation option to save a certain number of frames (png files) into a folder called 'Frames'. Execute the Python script `create_mp4.py`. The script will convert the frames from the 'Frames' folder into an MP4 animation. Once the script finishes running, you will find the `burst.mp4` video in the 'Output' folder

### Note on Data Accuracy

It is important to note that while the main geometry and input data are similar to the real experimental data, they are not exact replicas. For more accurate simulations, it is recommended to request the data directly from the NEA (Nuclear Energy Agency).