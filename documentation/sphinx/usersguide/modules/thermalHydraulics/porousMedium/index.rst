.. _userguide_thermalhydraulics_porousMedium:


The *phaseProperties* dictionary
--------------------------------

The *phaseProperties* dictionary can be found in *constant/fluidRegion/*. It is
a large dictionary that can be used to: choose the sub-solver to be used
(one-phase, legacy one-phase or two-phase); set various properties of the phases
(besides basic thermo-physical properties defined in the
*thermophysicalProperties* dictionary); set the properties of the sub-scale
structures (fuel pins, heat exchangers, etc) in the porous zones, including the
possibility to assign a *powerModel* for power production (e.g., nuclear fuel,
or constant power) and the *passiveProperties* of another sub-structure that
interacts thermally with the fluid (for instance the wrappers in sodium fast
reactors). In addition, models are available to model pumps and heat exchangers.
The name of the porous zones must coincide with that of the cellZones of the
fluidRegion mesh.

One can find detailed, commented examples in the tutorials
`3D_SmallESFR <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/3D_SmallESFR/extendedThermoMechanics/constant/fluidRegion/phaseProperties>`_ (single phase) and
`1D_boiling <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/featureCases/1D_boiling/constant/fluidRegion/phaseProperties>`_
(two phases). In addition, an example of how to use a two-dimensional
flow-regime map can be found in
`1D_PSBT_SC <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/featureCases/1D_PSBT_SC/Phase_Ex1_12223/constant/fluidRegion/phaseProperties>`_.





.. toctree::
   :maxdepth: 2

    FFdrag