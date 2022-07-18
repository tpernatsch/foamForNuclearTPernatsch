# User manual {#USERMAN}

**Work in progress!!**

<div class="border-box" style='padding:0.1em; margin-left: 4em;  margin-right: 8em;  border: 1px solid gray; background-color:#f2bfbb; color:#05134a'>
<b>Important note on how to use this manual</b>

This short user manual is meant to provide the user with enough information to understand the logic and various options of GeN-Foam. However, consistent with the Doxygen philosophy, and with the objective of minimizing inconsisitencies between documentation and cource code, detailed usage information of sub-solvers, behavioural models, etc., are included directly in the header (.H) file of the corresponsing classes. This manual normally provides links to these header files. Otherwise, they can be searched using the search function at the top right of the page. The links provided in the documentation will bring you to the Doxygen page of the corresponding header file. In this page, one can find the mentioned usage information under the section *Classes*. 
<br><br> In a similar fashion, to describe the use of complex dictionaries (i.e., input filed), this manual provides link to one or more commented dictionaries that are available in the tutorials.
</p>
</div>

## Some theoretical background

A fairly general theoretical presentation of GeN-Foam is provided in Ref. \cite FIORINA201524. It is recommended to go through this paper before starting to use GeN-Foam. However, the paper is getting quite old and it is recommended to refer to Ref. \cite FIORINA2016212  for the diffusion solver, Ref. \cite FIORINA2017419 for the SP3 solver, Refs. \cite Fiorina2019DetailedOpenFoam \cite Fiorina2015ApplicationCodes for the thermal-mechanic solver and its use for mesh deformation, Ref. \cite Fiorina2019DetailedOpenFoam for the SN solver, and Refs. \cite Radman2019ADesign \cite RADMAN2021111178 \cite RADMAN2021111422 for single- and two-phase thermal-hydraulics.

Here below a couple of essential points.

**The multi-region approach**

GeN-Foam employs a traditional OpenFOAM multi-region approach to model different physics using different meshes. This implies that the *0*, *constant* and *system* folders of each case contains multiple folders, one for each physics. In particular, the regions *fluidRegion*, *neutroRegion* and *thermoMechanicalRegion* are employed in GeN-Foam for thermal-hydraulics, neutronics and thermal-mechanics. Please notice that a dummy mesh must always be present for each physics, even if not solved for. One may use the dummy meshes and initial/boundary conditions provided in the EMPTY tutorial case.


**The multi-zone approach**

In order to assign different properties (for instance, different porous medium properties or different cross sections) to different zones in a mesh, GeN-Foam employs the OpenFOAM concept of cellZone. Each mesh should then be divided in different cellZones. Each cellZone is associated to a name and this name is used in *constant/fluidRegion/phaseProperties*, *constant/neutroRegion/nuclearData*, *constant/themoMechanicalRegion/themoMechanicalProperties* to associate each cellZone with a set of properties. The creation of cellZones is normally allowed by all meshers, though different names are normally used (for instance, *physical entities* in gmsh and *groups* in Salome). In same cases, conversion of the mesh into an OpenFOAM format creates cellSet instead of cellZones. In these case, one can use the topoSet utility to convert cellSets into cellZones.

## In practice

Most of the capabilities of GeN-Foam are encapsulated into  different C++ classes (or sub-libraries) for *neutronics*, *thermalHydraulics* and *thermoMechanics*.

In addition, a 4th class called *multiPhysicsControl* is employed to streamline the handling of various multi-physics options.

GeN-Foam itself is nothing but a complex fairly coupling loop that calls various functionalities of the *neutronics*, *thermalHydraulics*, and *thermoMechanics* classes and takes care of transfering coupling fields among them. 

The classes *neutronics*, *thermalHydraulics*, and *thermoMechanics*  translate for the user into the 3 different regions (and meshes). There is no requirement for the three meshes to occupy the same region of space. Consistent mapping of fields is performed and a reference value is given to a field if no correspondence is found in the mesh where its value is being projected from.

<div class="border-box" style='padding:0.1em; margin-left: 4em;  margin-right: 8em;  border: 1px solid gray; background-color:#f2f3fa; color:#05134a'>
<b>The meshes</b>

A dummy mesh must always be present in all physics (region) directories, even if not solved for. The EMPTY case is already provided with minimal dummy meshes and consistent fields in the “0” folder. Be careful! In case of parallel calculations all your meshes will have to have a number of cells equal or higher than the number of domains you are decomposing your geometry into. In case you need more cells than what available in the EMPTY case, you can run a refineMesh
</p>
</div>

The following sections describe the use of the 3 main classes of GeN-Foam, incuding some theory and useful references:
* [*neutronics*](@ref NEUTRONICS)
* [*thermalHydraulics*](@ref TH)
* [*thermoMechanics*](@ref TM)

The following section describes instead the coupling strategy and the general GeN-Foam options.

* [*Coupling options and time stepping*](@ref COUPLING)




© All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 2021