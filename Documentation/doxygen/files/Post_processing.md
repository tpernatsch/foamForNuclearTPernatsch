# Postprocessing {#POSTPROCESSING}

Postprocessing can be performed using Paraview, the standard post-processing tool used with OpenFOAM.  Paraview is launched using the command line:

`paraFoam -region (region name)`

where the region name is fluidRegion, neutroRegion or thermalMechanicalRegion (without the parenthesis!).
In case of parallel calculations, one should first reconstruct each one of the three meshes using the command

`recontructPar -region (region name)`

where the region name is once again fluidRegion, neutroRegion or thermalMechanicalRegion.

NB: Please notice that paraFoam is not always automatically installed together with OpenFOAM. In addition, depending on the operating system, installation of paraFoam may turn out to be problematic. In these cases, one may try to install OpenFOAM using docker. 

Beside paraFoam, GeN-Foam also creates, in the case folder (or in the processor folder for parallel simulations),  the GeN-Foam.dat file that summurizes few main quantity of interests: time(s), keff(-), power(W), flux0 (m-2s-1), TFuel_Max, TFuel_Avg TFuel_Min, TCladding_Max, TCladding_Avg, TCladding_Min.

In addition, OpenFOAM  allows to use [function objects](https://www.openfoam.com/documentation/guides/latest/doc/guide-function-objects.html) to extract specific information during or after the simulation.  Also in this case it is essential to indicate which region you want the function object to be applied to, for instance:

`postProcess -func singleGraph -region neutroRegion`

