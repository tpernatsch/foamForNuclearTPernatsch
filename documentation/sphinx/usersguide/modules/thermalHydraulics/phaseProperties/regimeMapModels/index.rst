.. _modules_thermalHydraulics_porousMedium_*regimeMapModels:

The *regimeMapModels* sub-dictionary
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Regime maps
-----------

In foamForNuclear, it is possible to employ 1- and 2-dimensional regime maps to use
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

Examples of



The oneParameter type is a 1-D regime map.
Here, the map name is lamTurb and it depends on the fluid Reynolds number.
Mismatches in the regime bounds of adjacent regimes automatically create
an interpolation region, wherein the interpolation type can be changed via
the interpolationMode keyword, which currently offers either linear or
quadratic interpolation. In the oneParameter map, the lowest bound of the
"lowest" regime is extended to -inf and the upper bound of the "upmost"
regime is extended to +inf

The twoParameter type is a 2-D regime map.

To use a regimeMap in a model, select the "byRegime" type for said model.
All physicsModels support it. The name of the regime map to be used
is specified via a regimeMap keyword. Then, the model for each regime is
specfied in subDicts named like each corresponding regime in the regimeMap.
To see an example of this, have a look at the heatTransferModel used for
this case.


.. code :: cpp



