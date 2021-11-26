# GeN-Foam Theory

A fairly general theoretical presentation of GeN-Foam is provided in Ref. \cite FIORINA201524. It is recommended to go through this paper before starting to use GeN-Foam. However, the paper is getting quite old and it is recommended to refer to Ref. \cite FIORINA2016212  for the diffusion solver, Ref. \cite FIORINA2017419 for the SP3 solver, Refs. \cite Fiorina2019DetailedOpenFoam \cite Fiorina2015ApplicationCodes for the thermal-mechanic solver and its use for mesh deformation, Ref. \cite Fiorina2019DetailedOpenFoam for the SN solver, and Refs. \cite Radman2019ADesign \cite RADMAN2021111178 \cite RADMAN2021111422 for single- and two-phase thermal-hydraulics.

Here below a few essential points.

**The multi-region approach**

GeN-Foam employs a traditional OpenFOAM multi-region approach to model different physics using different meshes. This implies that the *0*, *constant* and *system* folders of each case contains multiple folders, one for each physics. In particular, the regions *fluidRegion*, *neutroRegion* and *thermoMechanicalRegion* are employed in GeN-Foam for thermal-hydraulics, neutronics and thermal-mechanics. Please notice that a dummy mesh must always be present for each physics, even if not solved for. One may use the dummy meshes and initial/boundary conditions provided in the EMPTY tutorial case.


**The multi-zone approach**

In order to assign different properties (for instance, different porous medium properties or different cross sections) to different zones in a mesh, GeN-Foam employs the OpenFOAM concept of cellZone. Each mesh should then be divided in different cellZones. Each cellZone is associated to a name and this name is used in *constant/fluidRegion/phaseProperties*, *constant/neutroRegion/nuclearData*, *constant/themoMechanicalRegion/themoMechanicalProperties* to associate each cellZone with a set of properties. The creation of cellZones is normally allowed by all meshers, though different names are normally used (for instance, *physical entities* in gmsh and *groups* in Salome). In same cases, conversion of the mesh into an OpenFOAM format creates cellSet instead of cellZones. In these case, one can use the topoSet utility to convert cellSets into cellZones.


**Thermal-hydraulics**

Both single- and two-phase simulations can be performed using GeN-Foam. The selection is made by changing the keyword *thermalHydraulicType* in *constant/fluidRegion/phaseProperties*. Both solvers were developed for a coarse-mesh porous-medium treatment of complex structures such as core and heat exchanger, and for a standard RANS treatment of clear-fluid regions. A coarse-mesh porous-medium treatment of the core implies that the core is modeled without resolving each pin. As a matter of fact, in principle and for consistency, the finest radial mesh chosen by a user should not finer than one cell per pin cell. Pressure drops are modeled based on user-selectable correlations. Correlations are also  used to evaluate the heat exchange with the fuel. The fuel can be modeled using a 1-D sub-scale model that will calculate temperatures in the cladding and pellet for the pins that are contained in a cell. In addition to fuel, one can model simpler structures. In this case, one single temperature (and not a radial profile) is evaluated for each cell. This can be used to model structures like the assembly wrappers, the reflectors, or for a simplified modelling of a heat exchanger. 

The dictionary *constant/fluidRegion/phaseProperties* must include all information about pressure drop and heat transfer correlations, as well indication about the sub-scale models to be used for each mesh zone.

With regards to the boundary conditions, OpenFOAM provides most of the boundary conditions one may need for thermal-hydraulics models. In addition, a few boundary conditions have been included in GeN-Foam and can be found in *GeN-Foam/classes/thermalHydraulics/src/boundaryConditions*. In particular, information on the use of each boundary condition can be found in the header files (.H). 


**Neutronics**

GeN-Foam provides 5 different models for neutronics, namely: point kinetics, diffusion (transient and eigenvalue), adjoint diffusion (only eigenvalue), SP3 (transient and eigenvalue), discrete ordinates (only eigenvalue). One can select the desired model in  *constant/neutroRegion/neutronicsProperties*.

Nuclear data must be provided in *constant/neutroRegion/nuclearData*. The dictionaries *nuclearData...* give the possibility to parametrize the cross sections by giving a different set of cross sections in a perturbed state. GeN-Foam will perform a cell-by-cell linear interpolation except for fuel temperature, in which case a square root or logarithmic interpolation is employed
depending on the keyword *fastNeutrons* in *constant/neutroRegion/nuclearData*.

All GeN-Foam neutronics models can be used for liquid-fuel reactors, except for the point kinetics model. One can activate this option using the  *liquidFuel* keyword in */system/controlDict*. Of course, in such case one should pay attention to setting proper boundary conditions for the precursors.

Standard OpenFOAM boundary conditions such as *zeroGradient* and *fixedValue* can be employed. In addition the *albedoSP3* boundary condition is available in GeN-Foam and can be used both for the diffusion and SP3 solver. Details about the *albedoSP3* boundary condition can be found in its header file in *GeN-Foam/classes/neutronics/albedoSP3*. 

A specificity of the neutronic solvers is that GeN-Foam allows for the possibility of selecting the same boundary conditions for all neutron fluxes (in energy and angle), and the same boundary conditions for all the precursors groups. This can be done by indicating these boundary conditions in the *defaultFlux* and *defaultPrec* files. Of course, one may provide different boundary conditions for different fluxes and precursors by providing the related files in the */0/* folder (or in the folder related to the time step the simulation will start from). For example, in a 6-group diffusion simulation, one may set a general boundary for all energies using defaultFlux, and specific boundary conditions for the first and second group by: copying twice *defaultFlux* as *fluxStar0* and *fluxStar1*; and applying different boundary conditions to these two fluxes. 

NB: boundary conditions must be applied to *fluxStar...* and not to *flux...* since GeN-Foam solves for these variables. *fluxStar...* represent continuous fluxes, while *flux...* represent the real fluxes. They differ only in case discontinuity factors are employed /cite FIORINA2016212. 

NB2: defaultPrec has 1/m3 units except for the adjoint solver that needs 1/m2/s.


**Thermal-mechanics**
The thermal-mechanics solver of GeN-Foam is a simple linear elasticity solver that can be used to evaluate thermal deformations in a core. Temperatures are projected from the thermal-hydraulic solver (temperatures of fuels and structures) and the deformation field is used to deform the mesh for neutronics and thermal-hydraulics. In particular, the radial deformation of structures and the axial deformation of fuel are employed. to deform the neutronics mesh. 

As regards boundary conditions, besides the standard ones available in OpenFOAM, GeN-Foam includes a *tractionDisplacement* boundary condition that allows to set a pressure or a traction on a boundary. Work is ongoing to adapt to GeN-Foam the contact boundary condition of the OFFBEAT fuel behavior solver \cite SCOLARO2020110416. 

**Multi-physics coupling**
GeN-Foam can selectively couple each physics. The selection of the physics to solve is done in *system/controlDict* while the parameters for the coupling are set in *system/fvSolution*.










