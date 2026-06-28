.. _userguide_thermalhydraulics_turbulenceProperties:

The *turbulenceProperties* dictionary
-------------------------------------

*turbulenceProperties* is a standard OpenFOAM dictionary that allows you to define the turbulence model(s) to be used. For imported solvers such as :ref:`rhoPimpleFoam <rhoPimpleFoam>`, you can use the turbulence models that are usually available. For solvers that provide one- or two-phase porous-medium capabilities ( :ref:`onePhase <onePhase>`, :ref:`onePhaseLegacy <onePhaseLegacy>`, :ref:`twoPhase <twoPhase>`), foamForNuclear provides the models listed below.

When porous zones are present in a one-phase simulation, it is recommended to use *porousKEpsilon* (see :ref:`porousKEpsilon <porousKEpsilon>`). The Lahey model (see :ref:`LaheyKEpsilon <LaheyKEpsilon>`) and a mixture model (see :ref:`mixtureKEpsilon <mixtureKEpsilon>`) can be used for clear-fluid cases, or for mixed clear-fluid and porous-medium simulations in strongly advective two-phase flow scenarios where turbulent mixing may be neglected. In addition, a simple extension of the *porousKEpsilon* model has been implemented. This extension corrects the turbulent intensity using a term proportional to the fraction of the other phase (see :ref:`porousKEpsilon2PhaseCorrected <porousKEpsilon2PhaseCorrected>`).


.. raw:: html

   <br><br>


.. note ::
  Multiphase simulations require a ``turbulenceProperties`` dictionary for every phase (for example, ``turbulenceProperties.liquid`` and ``turbulenceProperties.vapour``).