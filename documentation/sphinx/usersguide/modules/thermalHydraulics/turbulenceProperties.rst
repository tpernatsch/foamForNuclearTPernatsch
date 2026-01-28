.. _userguide_thermalhydraulics_turbulenceProperties:

The *turbulenceProperties* dictionary
-------------------------------------

*turbulenceProperties*  is a standard OpenFOAM dictionary that allows defining the turbulence model(s) to be used. For imported solvers like :ref:`rhoPimpleFoam <rhoPimpleFoam>`, one can use the usually available turbulence models. For solvers that feature one- or two-phase porous-medium capabilities ():ref:`onePhase <onePhase>`, :ref:`onePhaseLegacy <onePhaseLegacy>`, :ref:`twoPhase <twoPhase>`), foamForNuclear provides the models reported below. When porous zones are present in a one-phase simulation, it is recommended to use
*porousKEpsilon* (see :ref:`porousKEpsilon.H <porousKEpsilon>`). The Lahey model (see
:ref:`LaheyKEpsilon.H <LaheyKEpsilon>`) and a mixture model (see
:ref:`mixtureKEpsilon.H <mixtureKEpsilon>`) can be used for
clear-fluids, or mixed clear-fluid and porous-medium simulations in case of
strongly advective two-phase flow scenarios where turbulent mixing may be
neglected. In addition, a simple extension of the *porousKEpsilon* model has
been implemented that allows to correct the turbulent intensity using a term
that is proportional to the fraction of the other phase (see
:ref:`porousKEpsilon2PhaseCorrected.H <porousKEpsilon2PhaseCorrected>`).

toctreeHere


.. note ::
  Multiphase simulation require a ``turbulenceProperties`` dictionary for every phase (e.g., turbulenceProperties.liquid and turbulenceProperties.vapour)