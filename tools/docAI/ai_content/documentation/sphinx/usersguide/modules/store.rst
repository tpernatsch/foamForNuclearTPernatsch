===============
Document Title
===============   
(Level 1)

---------------
Chapter Title
---------------   
(Level 2)

Section Title
=============   
(Level 3)

Subsection
----------   
(Level 4)

Sub-subsection
~~~~~~~~~~~~~   
(Level 5)




Physical properties (the ``constant`` folders)
==============================================

All the data for the GeN-Foam simulations can be provided in the following
input files (dictionaries):

- ``constant/thermoMechanicalRegion/thermoMechanicalProperties``: thermo-mechanical properties of structures, subdivided according to the ``cellZones`` of the ``thermoMechanicalRegion`` mesh. You can find a detailed, commented example in the tutorial `3D_SmallESFR <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/3D_SmallESFR/extendedThermoMechanics/constant/thermoMechanicalRegion/thermoMechanicalProperties>`_.
- ``constant/fluidRegion/g``: gravitational acceleration.
- ``constant/fluidRegion/turbulenceProperties``: standard OpenFOAM dictionary used to define the turbulence model. You can find a detailed, commented example in the tutorial `3D_SmallESFR <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/3D_SmallESFR/extendedThermoMechanics/constant/fluidRegion/turbulenceProperties>`_.
- ``constant/fluidRegion/thermophysicalProperties`` (for single-phase simulations): standard OpenFOAM dictionary used to define the thermo-physical properties of the coolant. You can find a detailed, commented example in the tutorial `3D_SmallESFR <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/3D_SmallESFR/extendedThermoMechanics/constant/fluidRegion/thermophysicalProperties>`_ (single phase).
- ``constant/fluidRegion/thermophysicalProperties.(name of fluid)`` (for two-phase simulations): standard OpenFOAM dictionaries used to define the thermo-physical properties of the different phases. The name of the fluid is defined in ``constant/fluidRegion/phaseProperties``. You can find a detailed, commented example in the tutorial `1D_boiling (liquid) <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/featureCases/1D_boiling/constant/fluidRegion/thermophysicalProperties.liquid>`_, `(vapour) <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/featureCases/1D_boiling/constant/fluidRegion/thermophysicalProperties.vapour>`_.
- ``constant/fluidRegion/phaseProperties``: a large dictionary that can be used to: determine whether the simulation is single-phase or two-phase; set various properties of the phases (in addition to the thermo-physical properties defined in ``constant/fluidRegion/thermophysicalProperties``); set the properties of the sub-scale structures (fuel pins, heat exchangers, etc.) in the porous zones, including the possibility to assign a ``powerModel`` for power production (e.g., nuclear fuel, or constant power) and the ``passiveProperties`` of another sub-structure that interacts thermally with the fluid (for instance, the wrappers in sodium fast reactors). The name of the porous zones must match that of the ``cellZones`` of the ``fluidRegion`` mesh. Anisotropic pressure drops can be set using the keywords ``transverseDragModel`` (Blasius, GunterShaw, same) and ``principalAxis`` (localX, localY, localZ) in the sub-dictionary ``dragModels.(nameOfPhase).structure.(nameOfCellZones)``. ``principalAxis`` sets the axis on which the nominal ``dragModel`` is used. ``transverseDragModel`` sets the model to be used in the two directions perpendicular to ``principalAxis``. If ``same`` is chosen as ``transverseDragModel``, the code uses the nominal model in all directions, while allowing an anisotropic hydraulic diameter. The anisotropy of the hydraulic diameter can be set using the keyword ``localDhAnisotropy`` and assigning it a vector of 3 scaling factors (one for each local direction). You can find detailed, commented examples in the tutorials `3D_SmallESFR <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/3D_SmallESFR/extendedThermoMechanics/constant/fluidRegion/phaseProperties>`_ (single phase) and `1D_boiling <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/featureCases/1D_boiling/constant/fluidRegion/phaseProperties>`_ (two phases).
- ``constant/neutroRegion/neutronicsProperties``: dictionary used to control whether the calculation is an eigenvalue calculation or a transient. You can find detailed, commented examples in most tutorials. See, for instance, `3D_SmallESFR <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/3D_SmallESFR/extendedThermoMechanics/constant/neutroRegion/neutronicsProperties>`_ (single phase).
- ``timeStep/uniform/reactorState``: contains the target power (``pTarget``) for eigenvalue calculations, the ``keff`` resulting from the eigenvalue calculations, and the external reactivity (i.e., the additional reactivity that can be added, for instance, to simulate a reactivity step). N.B.: ``keff`` has no effect on pointKinetics. You can find detailed, commented examples in most tutorials. N.B.2: In point kinetics, ``pTarget`` is the initial value used by the point kinetics solver to plot results. However, the solver actually scales the ``powerDensity`` and ``flux`` fields provided by the user. It is up to the user to ensure that ``pTarget`` is consistent with the provided ``powerDensity`` and ``flux`` fields. A commented ``reactorState`` can be found in `3D_SmallESFR <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/3D_SmallESFR/extendedThermoMechanics/constant/neutroRegion/reactorState>`_ (single phase). Please note that eigenvalue calculations update the ``keff`` value in this dictionary. In parallel calculations, the updated value can be found in ``processor0/constant/neutroRegion/reactorState``.
- ``constant/neutroRegion/nuclearData``: contains all basic nuclear properties for the reference and perturbed reactor states. In addition, this file includes information about the perturbed and reference parameters. For instance, perturbing the fuel temperature must include ``TFuel`` at the reference state and any number of perturbed temperatures. GeN-Foam performs interpolation between the reference and perturbed reactor states. If no data are provided, the reference cross-sections are used. Nuclear data can be generated using any nuclear code. Several tools are provided with GeN-Foam (in the ``tools`` folder) that automatically convert Monte Carlo output files into the nuclear data files used by GeN-Foam (see ``tools``). The entry ``discFactor`` is used only if discontinuity factors must be applied. The term ``integralFlux`` is used only if the automatic adjustment of discontinuity factors is performed (see :ref:`FIORINA2016212 <FIORINA2016212>`). Nonetheless, these entries should always be present. More information is available in :ref:`nuclearData section <userguide_neutronics_nuclearData>`. You can find detailed, commented examples of ``nuclearData`` in the tutorials `3D_SmallESFR <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/3D_SmallESFR/extendedThermoMechanics/constant/neutroRegion/nuclearData>`_ (for diffusion or SP3), `Godiva_SN <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/Godiva_SN/constant/neutroRegion/nuclearData>`_ (for discrete ordinates), and `2D_onePhaseAndPointKineticsCoupling <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/featureCases/2D_onePhaseAndPointKineticsCoupling/rootCase/constant/neutroRegion/nuclearData>`_ (for point kinetics). You can also find examples of the ``nuclearData...`` files in the tutorial `3D_SmallESFR <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/3D_SmallESFR/extendedThermoMechanics/constant/neutroRegion>`_.
- ``constant/neutroRegion/quadratureSet``: contains the quadrature set for discrete ordinate calculations (see :ref:`Quadrature set dictionary <userguide_neutronics_quadratureSet>` for more details). You can find examples of three different quadrature sets in the tutorial `Godiva_SN <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/Godiva_SN/constant/neutroRegion/>`_.
- ``constant/neutroRegion/CRMove``: contains input data for control rod movement. Control rods can be moved from the initial position to a new one by selecting an initial and final time for the insertion/extraction, and the insertion/extraction speed (positive speed for insertion) (see :ref:`Control rod movement <userguide_neutronics_crmove>` for more details). You can find a commented example in the tutorial `3D_SmallESFR <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/3D_SmallESFR/extendedThermoMechanics/constant/neutroRegion/CRmove>`_, although this option is not actually used in the tutorial.




Initial values and boundary conditions (the time folders)
=========================================================

As in all standard OpenFOAM solvers, initial values (IC) and boundary
conditions (BC) should be provided in the ``0`` folder, or in the time folder
corresponding to the ``startTime`` of the simulation, if it is different from
``0``.

Fields that need to be mapped across physics can be initialized by running
``GeN-Foam -initializeMappedFields``. This command constructs the solvers
and mappings, and maps the fields across physics as specified in the
``regionsDict``. This might be useful for restart simulations and for
point-kinetics simulations, as explained in
:ref:`this page <userguide_thermalhydraulics_settingcase_setpower>`.




Discretization and solution (the ``system`` folders) 
====================================================

Details for the discretization and solution of single-physics equations are
handled in a standard OpenFOAM way, i.e., through the ``fvSolution`` and
``fvSchemes`` dictionaries in ``system/regionName`` (e.g., ``system/fluidRegion``).
Coupling and simulation details are determined through the ``controlDict`` in
the ``system`` folder. The ``controlDict`` is significantly extended compared
to a standard OpenFOAM ``controlDict`` in order to control which solvers are
activated and the flags that affect the behavior of GeN-Foam as a whole.