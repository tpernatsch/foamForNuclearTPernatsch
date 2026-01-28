.. _userguide_thermalhydraulics_turbulenceProperties:

The *turbulenceProperties* dictionary
-------------------------------------


*turbulenceProperties*  is a standard OpenFOAM dictionary that allows defining the turbulence model to be used. In addition to the several models already available in OpenFOAM, foamForNuclear provides some specialized models.


When porous zones are present in the simulation, it is recommended to use
*porousKEpsilon* (see :ref:`porousKEpsilon.H <porousKEpsilon>`). The only difference w.r.t. the
standard k-epsilon model is that it forces k and epsilon to equilibrium values
inside the porous zones. These equilibrium values can be set in the
*porousKepsilonProperties* sub-dictionary. Please notice that a porous medium
simulation using the equilibrium values of k and epsilon for the sub-scale
structure (viz., the values inside a fuel sub-channel) would entail the risk of
an unstable solution. This occurs because the turbulent viscosity is primarily
associated with the sub-scale structure, which might not be sufficient to
maintain stability at the coarse mesh's length scale. To address this problem,
one can define the keyword DhStruct in
*constant/fluidRegion/phaseProperties/dragModels.(nameOfPhase).structure.(nameOfCellZones)*.
This keyword defines the hydraulic diameter of the whole porous structure (viz.,
the dimension of the assembly, if using baffles to model wrappers, or of the
entire core). The code uses it to make sure the turbulent viscosity results in a
laminar Reynolds number (defaulted to 500).

While some approaches to model k and epsilon for two-phase flow simulations are
presently included in the code. In particular, the Lahey model (see
:ref:`LaheyKEpsilon.H <LaheyKEpsilon>`) and a mixture model (see
:ref:`mixtureKEpsilon.H <mixtureKEpsilon>`) can be used for
clear-fluids, or mixed clear-fluid and porous-medium simulations in case of
strongly advective two-phase flow scenarios where turbulent mixing may be
neglected. In addition, a simple extension of the *porousKEpsilon* model has
been implemented that allows to correct the turbulent intensity using a term
that is proportional to the fraction of the other phase (see
:ref:`porousKEpsilon2PhaseCorrected.H <porousKEpsilon2PhaseCorrected>`).

One can find a detailed, commented example of a porous one-phase simulation in
the tutorial
`3D_SmallESFR <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/3D_SmallESFR/extendedThermoMechanics/constant/fluidRegion/turbulenceProperties>`_.

