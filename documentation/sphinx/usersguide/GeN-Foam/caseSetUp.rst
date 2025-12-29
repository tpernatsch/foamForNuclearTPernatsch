=============
Preprocessing
=============

Before running GeN-Foam, one has to provide meshes, physical properties,
discretization methods (if one does not want to use the default ones) and
simulation details. This section provides a quick overview of the input deck of
GeN-Foam.


Meshing
=======

GeN-Foam uses any number of meshes for that can be e.g neutronics,
thermal-hydraulics and thermal-mechanics. There is no requirement for the meshes
to occupy the same region of space. Consistent mapping of fields is performed
and a reference value is given to a field if no correspondence is found in the
mesh where its value is being projected from. It follows that the geometry for
neutronics can cover only a small part of the overall reactor geometry. Meshes
can be created with every OpenFOAM-compatible tool (e.g blockMesh, ANSYS/FLUENT,
...). Meshes can (should) be divided into zones (cellZones) to allow the use of
different physical properties (e.g., cross-sections, power models, ...) in
different reactor regions. Sometimes, when converting a mesh to the OpenFOAM
(``polyMesh`` folder) format, cellSets (and not cellZones) are created. The
``topoSetDict`` can be used to convert cellSets to cellZones.

The `3D_SmallESFR <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/3D_SmallESFR/>`_
tutorial includes an example of mesh generation with Gmsh `GMSH2009 <https://onlinelibrary.wiley.com/doi/10.1002/nme.2579>`_.
The folder includes the subfolder *3dESFRMesh* containing three subfolders for the
generation of the meshes for neutronics, thermal-hydraulics and
thermal-mechanics. To create the 3D ESFR core geometry, one should execute the
following command in a terminal:

.. code :: bash

    gmsh esfMain.geo


This will create (after a relatively long time) a geometry and open the Gmsh
graphical interface. The geometry must then be meshed and the resulting mesh
saved and copied in the root of the *3D_SmallESFR* folder.

By typing in a terminal:

.. code :: bash

    gmshToFoam nameMeshFile


A *polyMesh* folder will be created (or updated) in the folder *constant*. One
should then copy this folder in *constant/neutroRegion*, *constant/fluidRegion*
or *constant/thermoMechanicalRegion*, and repeat the operation for all meshes.

Please notice that the *3D_SmallESFR* tutorial already contains the correct
*polyMesh* folders so that one can avoid the mesh generation step.


Physical properties
===================

All the data for the GeN-Foam simulations can be filled in the following input
files (dictionaries):

- ``constant/thermoMechanicalRegion/thermoMechanicalProperties``: thermo-mechanical properties of structures, subdivided according to the cellZones of the thermoMechanicalRegion mesh. One can find a detailed, commented example in the tutorial `3D_SmallESFR <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/3D_SmallESFR/extendedThermoMechanics/constant/thermoMechanicalRegion/thermoMechanicalProperties>`_.
- ``constant/fluidRegion/g``: gravitational acceleration.
- ``constant/fluidRegion/turbulenceProperties``: standard OpenFOAM dictionary to define the turbulence model to be used. One can find a detailed, commented example in the tutorial `3D_SmallESFR <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/3D_SmallESFR/extendedThermoMechanics/constant/fluidRegion/turbulenceProperties>`_.
- ``constant/fluidRegion/thermophysicalProperties`` (for single-phase simulations): standard OpenFOAM dictionary to define the thermo-physical properties of the coolant. One can find a detailed, commented example in tutorial `3D_SmallESFR <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/3D_SmallESFR/extendedThermoMechanics/constant/fluidRegion/thermophysicalProperties>`_ (single phase)
- ``constant/fluidRegion/thermophysicalProperties.(name of fluid)`` (for two-phase simulations): standard OpenFOAM dictionaries to define the thermo-physical properties of various phases. The name of fluid is defined in *constant/fluidRegion/phaseProperties*. One can find a detailed, commented example in the tutorial `1D_boiling (liquid) <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/featureCases/1D_boiling/constant/fluidRegion/thermophysicalProperties.liquid>`_, `(vapour) <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/featureCases/1D_boiling/constant/fluidRegion/thermophysicalProperties.vapour>`_.
- ``constant/fluidRegion/phaseProperties``: a large dictionary that can be used to: determine whether the simulation is single-phase or two-phase; set various properties of the phases (besides the thermo-physical properties defined in *constant/fluidRegion/thermophysicalProperties*); set the properties of the sub-scale structures (fuel pins, heat exchangers, etc) in the porous zones, including the possibility to assign a *powerModel* for power production (e.g., nuclear fuel, or constant power) and the *passiveProperties* of another sub-structure that interacts thermally with the fluid (for instance the wrappers in sodium fast reactors). The name of the porous zones must coincide with that of the cellZones of the fluidRegion mesh. Anisotropic pressure drops can be set by using the keywords *transverseDragModel* (Blasius, GunterShaw, same) and *principalAxis*(localX, localY, localZ) in the sub-dictionary *dragModels.(nameOfPhase).structure.(nameOfCellZones)*. *principalAxis* sets the axis on which the nominal dragModel is used. *transverseDragModel* sets the model to be used in the two directions that are perpendicular to *principalAxis*. If *same* is chosen as *transverseDragModel*, the code will use the nominal model in all directions, but with the possibility of an anisotropic hydraulic diameter. The anisotropy of the hydraulic diameter can be set using the keyword *localDhAnisotropy* and assigned to it a vector of 3 scaling factors (one for each local direction). One can find detailed, commented examples in the tutorials `3D_SmallESFR <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/3D_SmallESFR/extendedThermoMechanics/constant/fluidRegion/phaseProperties>`_ (single phase) and `1D_boiling <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/featureCases/1D_boiling/constant/fluidRegion/phaseProperties>`_ (two phases).
- ``constant/neutroRegion/neutronicsProperties``: dictionary to control if it's an eigenvalue calculation or a transient. One can find detailed, commented examples in most tutorials. See for instance `3D_SmallESFR <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/3D_SmallESFR/extendedThermoMechanics/constant/neutroRegion/neutronicsProperties>`_ (single phase).
- ``timeStep/uniform/reactorState``: contains the target power (pTarget) for eigenvalue calculations, the keff that results from the eigenvalue calculations and the external reactivity (i.e., the extra reactivity one can add for instance to simulate a reactivity step). N.B.: keff has no effect on pointKinetics. You can find detailed, commented examples in most tutorials. N.B.2: In point kinetics, pTarget is the initial value used by the point kinetics solver to plot results, but the solver actually scales the *powerDensity* and *flux* fields provided by the user. It is up to the user to make sure that pTarget is consistent with the powerDensity and flux fields. A commented reactorState can be found in `3D_SmallESFR <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/3D_SmallESFR/extendedThermoMechanics/constant/neutroRegion/reactorState>`_ (single phase). Please note that eigenvalue calculations will update the keff value in this dictionary. In parallel calculations, the updated value can be found in *processor0/constant/neutroRegion/reactorState*.
- ``constant/neutroRegion/nuclearData``: contains all basic nuclear properties for the reference and perturbed reactor states. In addition, this file include information about the perturbed and reference parameters. For instance, perturbing the fuel temperature must include *TFuel* at reference and any number of perturbed temperatures. Interpolation is performed by GeN-Foam between reference and perturbed reactor states. If no data are provided, the reference cross-sections are used. Nuclear data can be generated using any nuclear code. Several tools are provided with GeN-Foam (in the *tools* folder) that automatically converts Monte Carlo output files into the nuclear data files employed by GeN-Foam (see *tools*). The entry *discFactor* is used only if discontinuity factors have to be used. The term *integralFlux*, is used only if the automatic adjustment of discontinuity factors is performed (see :ref:`FIORINA2016212 <FIORINA2016212>`). Nonetheless, these entries should always be present. More information in :ref:`nuclearData section <userguide_neutronics_nuclearData>`. One can find detailed, commented examples of nuclearData in the tutorials `3D_SmallESFR <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/3D_SmallESFR/extendedThermoMechanics/constant/neutroRegion/nuclearData>`_ (for diffusion or SP3), `Godiva_SN <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/Godiva_SN/constant/neutroRegion/nuclearData>`_ (for discrete ordinates) and `2D_onePhaseAndPointKineticsCoupling <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/featureCases/2D_onePhaseAndPointKineticsCoupling/rootCase/constant/neutroRegion/nuclearData>`_ (for point kinetics). One can find examples of the *nuclearData...* files in the tutorial `3D_SmallESFR <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/3D_SmallESFR/extendedThermoMechanics/constant/neutroRegion>`_
- ``constant/neutroRegion/quadratureSet``: contains the quadrature set for discrete ordinate calculations (see :ref:`Quadrature set dictionary <userguide_neutronics_quadratureSet>` for more details). One can find examples of three different quadrature sets in the tutorial `Godiva_SN <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/Godiva_SN/constant/neutroRegion/>`_.
- ``constant/neutroRegion/CRMove``: contains input data for control rods movement. Control rods can be moved from the initial position to a new one by selecting an initial and final time of the insertion/extraction and the speed of insertion/extraction (positive speed for insertion) (see :ref:`Control rod movement <userguide_neutronics_crmove>` for more details). One can find a commented example in the tutorial `3D_SmallESFR <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/3D_SmallESFR/extendedThermoMechanics/constant/neutroRegion/CRmove>`_, though this option is not actually used in the tutorial.


Initial values and boundary conditions
======================================

As in all standard OpenFOAM solvers, initial values (IC) and boundary conditions
(BC) should be provided in the “0” folder, or in the time folder corresponding
to the *startTime* of the simulation, if different than 0.

Fields which need to be mapped across physics can be initialized by running
``GeN-Foam -initializeMappedFields``. This will construct the solvers and
mappings and map the fields across physics as specified in the
*regionsDict*. This might be useful for restart simulations and for
point-kinetics simulations as explained in
:ref:`this page <userguide_thermalhydraulics_settingcase_setpower>`.


Discretization and solution
===========================

Details for discretization and solution of single-physics equations are handled
in a standard OpenFOAM way, i.e., through the *fvSolution* and *fvSchemes*
dictionaries in *system/regionName* (e.g *system/fluidRegion*). Coupling and
simulation details are determined through the *controlDict* in the *system*
folder. The *controlDict* is significantly extended compared to a standard
OpenFOAM controlDict in order to control what solvers are activated and flags
that affect the behavior of GeN-Foam as a whole.
