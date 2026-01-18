.. _modules_thermalHydraulics_porousMedium_regimeMapModels:

The *regimeMapModels* sub-dictionary
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is possible to employ 1- and 2-dimensional regime maps to use
different models for different flow conditions. This can be used for instance
in one-phase simulation to provide different correlations for turbulent and
laminar flow (see `3D_SmallESFR <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/3D_SmallESFR/extendedThermoMechanics/constant/fluidRegion/phaseProperties>`_
for a commented example), or in 2-phase flow simulations to provide full regime
maps (see `1D_boiling <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/featureCases/1D_boiling/constant/fluidRegion/phaseProperties>`_
for a commented example of a 1-dimensional map, and
`1D_PSBT_SC <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/featureCases/1D_PSBT_SC/Phase_Ex1_12223/constant/fluidRegion/phaseProperties>`_
for a commented example of a 2-dimensional map). Multiple maps can be used in
the same simulation. In addition, regime-maps models can be mixed with
multi-regime models, as in
`1D_PSBT_SC <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/featureCases/1D_PSBT_SC/Phase_Ex1_12223/constant/fluidRegion/phaseProperties>`_,
where a *preCHFTraceRegimeMap* is employed to assign models for phase
dispersion, interfacial area and bubble diameter, while a single multi-regime
model is employed to describe heat transfer between liquid and structure
throughout the various regimes.

An example of 1-D model is reported below.
Here, the map name is *lamTurb* (an arbitray name) and the type *oneParameter*. To use a regimeMap in a model, the user should select the "byRegime" type for said model. All :ref:`physics models <modules_thermalHydraulics_porousMedium_physicsModels>` support it. The name of the regime map to be used
is specified via a regimeMap keyword. Then, the model for each regime is
specfied in subDicts named like each corresponding regime in the regimeMap. See :ref:`physics models <modules_thermalHydraulics_porousMedium_physicsModels>` for additional details.

.. code :: cpp
   regimeMapModels
   {
      "lamTurb"
      {
         type                oneParameter;
         ...
      }
   }

Avilable regime maps types include:













