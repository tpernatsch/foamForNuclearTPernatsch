.. _userguide_thermalhydraulics_thermophysicalProperties:

The *thermophysicalProperties* dictionary
-----------------------------------------

The *thermophysicalProperties* dictionary is located under
*constant/fluidRegion/*. It is a standard OpenFOAM dictionary that allows
you to define the thermophysical properties of the coolant. For two-phase
flow analyses, you must use two dictionaries named
*thermophysicalProperties.(name of fluid)*. The names of the two fluids are
defined in the *phaseProperties* dictionary.

Detailed examples are available in the tutorials:
`3D_SmallESFR <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/3D_SmallESFR/extendedThermoMechanics/constant/fluidRegion/thermophysicalProperties>`_ (one-phase),
`1D_boiling (liquid) <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/featureCases/1D_boiling/constant/fluidRegion/thermophysicalProperties.liquid>`_ (two-phase, liquid),
`1D_boiling (vapour) <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/featureCases/1D_boiling/constant/fluidRegion/thermophysicalProperties.vapour>`_ (two-phase, vapour).