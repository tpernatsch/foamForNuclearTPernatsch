.. _userguide_thermalhydraulics_thermophysicalProperties:

The *thermophysicalProperties* dictionary
-----------------------------------------

The *thermophysicalProperties* dictionary can be found under
*constant/fluidRegion/*. It is a standard OpenFOAM dictionary that allows
defining the thermo-physical properties of the coolant. When performing
two-phase flow analyses, two dictionaries must be employed named
*thermophysicalProperties.(name of fluid)*. The name of the two fluids is
defined in the *phaseProperties* dictionary.

One can find a detailed, commented example in the tutorials
`3D_SmallESFR <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/3D_SmallESFR/extendedThermoMechanics/constant/fluidRegion/thermophysicalProperties>`_ (one-phase),
`1D_boiling (liquid) <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/featureCases/1D_boiling/constant/fluidRegion/thermophysicalProperties.liquid>`_ (two-phase, liquid)
`1D_boiling (vapour) <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/featureCases/1D_boiling/constant/fluidRegion/thermophysicalProperties.vapour>`_ (two-phase, vapour)