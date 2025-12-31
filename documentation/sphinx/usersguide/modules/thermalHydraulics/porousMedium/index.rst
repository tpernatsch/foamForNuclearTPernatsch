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





.. list-table:: Fluid-structure drag models (:ref:`FSDragCoefficientModel.H <FSDragCoefficientModel>`)
    :widths: 50 50 50

    * - Baxi Dalle Donne
      - :ref:`BaxiDalleDonneFSDragCoefficient.H <BaxiDalleDonneFSDragCoefficient>`
      - N/A
    * - Churchill
      - :ref:`ChurchillFSDragCoefficient.H <ChurchillFSDragCoefficient>`
      - N/A
    * - Engel as in TRACE
      - :ref:`EngelFSDragCoefficient.H <EngelFSDragCoefficient>`
      - N/A
    * - Modified Engel
      - :ref:`modifiedEngelFSDragCoefficient.H <modifiedEngelFSDragCoefficient>`
      - N/A
    * - No Kazimi
      - :ref:`NoKazimiFSDragCoefficient.H <NoKazimiFSDragCoefficient>`
      - N/A
    * - Rehme
      - :ref:`RehmeFSDragCoefficient.H <RehmeFSDragCoefficient>`
      - N/A
    * - Drag coefficient as Reynolds power
      - :ref:`ReynoldsPowerFSDragCoefficient.H <ReynoldsPowerFSDragCoefficient>`
      - :math:`A \times Re^B+C`
    * - Drag coefficient as Colebrook correlation
      - :ref:`ColebrookFSDragCoefficient.H <ColebrookFSDragCoefficient>`
      - :math:`(A \log_{10}(Re)+B)^C`


.. list-table:: Two-phase drag multipliers (:ref:`twoPhaseDragMultiplierModel.H <twoPhaseDragMultiplierModel>`)
    :widths: 50 50

    * - Chen Kalish
      - :ref:`ChenKalishTwoPhaseDragMultiplier.H <ChenKalishTwoPhaseDragMultiplier>`
    * - Constant
      - :ref:`constantTwoPhaseDragMultiplier.H <constantTwoPhaseDragMultiplier>`
    * - Kaiser 74
      - :ref:`Kaiser74TwoPhaseDragMultiplier.H <Kaiser74TwoPhaseDragMultiplier>`
    * - Kaiser 88
      - :ref:`Kaiser88TwoPhaseDragMultiplier.H <Kaiser88TwoPhaseDragMultiplier>`
    * - Kottowski Savatteri
      - :ref:`KottowskiSavatteriTwoPhaseDragMultiplier.H <KottowskiSavatteriTwoPhaseDragMultiplier>`
    * - Lockhart Martinelli
      - :ref:`LockhartMartinelli.H <LockhartMartinelli>`
    * - Lottes Flinn
      - :ref:`LottesFlinnTwoPhaseDragMultiplier.H <LottesFlinnTwoPhaseDragMultiplier>`
    * - Lottes Flinn Nguyen
      - :ref:`LottesFlinnNguyenTwoPhaseDragMultiplier.H <LottesFlinnNguyenTwoPhaseDragMultiplier>`


Heat transfer models
--------------------

Currently, available models to describe the energy exchange with a sub-scale
structure or with a second phase include:

.. list-table:: Fluid-fluid heat-transfer models (:ref:`FFHeatTransferCoefficientModel.H <FFHeatTransferCoefficientModel>`)
    :widths: 50 50 50

    * - No Kazimi
      - :ref:`NoKazimiFFHeatTransferCoefficient.H <NoKazimiFFHeatTransferCoefficient>`
      - N/A
    * - Nusselt number correlation as Reynolds and Prandtl powers
      - :ref:`NusseltFFHeatTransferCoefficient.H <NusseltFFHeatTransferCoefficient>`
      - :math:`Nu = A + B \times Re^C Pr^D`


.. list-table:: Fluid-structure heat-transfer models (see FSHeatTransferCoefficientModel.H)
    :widths: 50 50 50

    * - Nusselt number correlation as Reynolds and Prandtl powers and surface to bulk temperature ratio
      - :ref:`NusseltFSHeatTransferCoefficient.H <NusseltFSHeatTransferCoefficient>`
      - :math:`Nu = A + B \times Re^C Pr^D \left( \frac{T_w}{T_b} \right)^E`
    * - Nusselt number correlation as Reynolds and Prandtl powers, plus an additional heat transfer coefficient to take into account the resistance of a wall
      - :ref:`NusseltAndWallFSHeatTransferCoefficient.H <NusseltAndWallFSHeatTransferCoefficient>`
      - :math:`H = \frac{Nu \times \kappa}{D_h} + H_{wall} \quad \text{with} \quad Nu = A + B \times Re^C Pr^D`
    * - Shah
      - :ref:`ShahFSHeatTransferCoefficient.H <ShahFSHeatTransferCoefficient>`
      - N/A
    * - Gorenflo
      - :ref:`GorenfloFSHeatTransferCoefficient.H <GorenfloFSHeatTransferCoefficient>`
      - N/A
    * - multiRegimeBoilingTRACE - multi-regime heat transfer coefficient that replicates what TRACE does below CHF (ADD REF GAUTHIER)
      - :ref:`multiRegimeBoilingTRACEFSHeatTransferCoefficient.H <multiRegimeBoilingTRACEFSHeatTransferCoefficient>`
      - N/A
    * - multiRegimeBoilingTRACE - multi-regime heat transfer coefficient that replicates what TRACE does, including CHF and post-CHF. Not verified! It requires specifying the *multiRegimeBoilingTRACECHF* model for the water.structure heat transfer, and the *multiRegimeBoilingVapourTRACE* model for the vapour.structure heat transfer. Please notice that a lookup table for CHF is still missing.
      - :ref:`multiRegimeBoilingTRACECHFFSHeatTransferCoefficient.H <multiRegimeBoilingTRACECHFFSHeatTransferCoefficient>` and :ref:`multiRegimeBoilingVapourTRACEFSHeatTransferCoefficient.H <multiRegimeBoilingVapourTRACEFSHeatTransferCoefficient>`
      - N/A
    * - multiRegimeBoiling - multi-regime heat transfer coefficient that replicates the same logic as TRACE, but with more flexibility for user-selectable sub-models. Can be used below CHF.
      - :ref:`multiRegimeBoilingFSHeatTransferCoefficient.H <multiRegimeBoilingFSHeatTransferCoefficient>`
      - N/A

- Sub-models employed by the multi-regime models:
    - Critical heat flux related models (:ref:`CHFModel.H <CHFModel>`):
        - Critical heat flux models
            - Constant, user-selectable value (:ref:`constantCHF.H <constantCHF>`)
            - Lookup table, not yet implemented (empty class at :ref:`lookUpTableCHF.H <lookUpTableCHF>`)
        - Leidenfrost models (:ref:`TLFModel.H <TLFModel>`)
            - Groeneveld Stewart (:ref:`GroeneveldStewartTLF.H <GroeneveldStewartTLF>`)
    - Flow Enhancement Factor Models (:ref:`flowEnhancementFactorModel.H <flowEnhancementFactorModel>`)
        - Chen (:ref:`ChenFlowEnhancementFactor.H <ChenFlowEnhancementFactor>`)
        - COBRA-TF (:ref:`COBRA-TFFlowEnhancementFactor.H <COBRA-TFFlowEnhancementFactor>`)
        - Rezkallah Sims (:ref:`RezkallahSimsFlowEnhancementFactor.H <RezkallahSimsFlowEnhancementFactor>`)
    - Post-CHF models
        - Cachard (for liquid) (:ref:`CachardLiquidFSHeatTransferCoefficient.H <CachardLiquidFSHeatTransferCoefficient>`)
        - Cachard (for vapour) (:ref:`CachardVapourFSHeatTransferCoefficient.H <CachardVapourFSHeatTransferCoefficient>`)
    - Sub-Cooled Boiling Fraction Models (:ref:`subCooledBoilingFractionModel.H <subCooledBoilingFractionModel>`)
        - Constant (:ref:`constantSubCooledBoilingFraction.H <constantSubCooledBoilingFraction>`)
        - Saha Zuber (:ref:`SahaZuberSubCooledBoilingFraction.H <SahaZuberSubCooledBoilingFraction>`)
    - Superposition Nucleate Boiling (:ref:`superpositionNucleateBoilingFSHeatTransferCoefficient.H <superpositionNucleateBoilingFSHeatTransferCoefficient>`)
    - Suppression factor models (:ref:`suppressionFactorModel.H <suppressionFactorModel>`)
        - Chen (:ref:`ChenSuppressionFactor.H <ChenSuppressionFactor>`)
        - COBRA-TF (:ref:`COBRA-TFSuppressionFactor.H <COBRA-TFSuppressionFactor>`)
    - Temperature of the onset of nucleate boiling (:ref:`TONBModel.H <TONBModel>`)
        - Basu (:ref:`BasuTONB.H <BasuTONB>`)


Special models
--------------

Dedicated models for specific sub-scale structures are described in more details
in :ref:`this page <thermalHydraulicsSubScaleStructures>`.


Models for two-phase flows
--------------------------

Currently, available models for two-phase flow simulations include:

- Contact partition models (:ref:`contactPartitionModel.H <contactPartitionModel>`)
    - Linear (:ref:`linearContactPartition.H <linearContactPartition>`)
    - Complementary (:ref:`complementaryContactPartition.H <complementaryContactPartition>`)
- Dispersions models (:ref:`dispersionModel.H <dispersionModel>`)
    - Constant (:ref:`constantDispersion.H <constantDispersion>`)
- Fluid diameter models (:ref:`fluidDiameterModel.H <fluidDiameterModel>`)
    - Iso-molar bubble (:ref:`isomolarBubbleFluidDiameter.H <isomolarBubbleFluidDiameter>`)
    - Iso-thermal bubble  (:ref:`isothermalBubbleFluidDiameter.H <isothermalBubbleFluidDiameter>`)
    - Pipe film (:ref:`pipeFilmFluidDiameter.H <pipeFilmFluidDiameter>`)
    - Wallis film (:ref:`WallisFilmFluidDiameter.H <WallisFilmFluidDiameter>`)
- Virtual mass models (:ref:`virtualMass.H <virtualMass>`)
    - Virtual mass coefficient (:ref:`virtualMassCoefficientModel.H <virtualMassCoefficientModel>`)
- Interfacial area models (:ref:`interfacialAreaModel.H <interfacialAreaModel>`)
    - Annular (:ref:`annularInterfacialArea.H <annularInterfacialArea>`)
    - No Kazimi (:ref:`NoKazimiInterfacialArea.H <NoKazimiInterfacialArea>`)
    - Schor (:ref:`SchorInterfacialArea.H <SchorInterfacialArea>`)
    - Spherical (:ref:`sphericalInterfacialArea.H <sphericalInterfacialArea>`)
- Phase change models
    - Forced constant (:ref:`forcedConstantPhaseChange.H <forcedConstantPhaseChange>`)
    - Heat driven (:ref:`heatDrivenPhaseChange.H <heatDrivenPhaseChange>`)
    - Latent heat (:ref:`latentHeatModel.H <latentHeatModel>`)
        - Fink Leibowitz for sodium (:ref:`FinkLeibowitzLatentHeat.H <FinkLeibowitzLatentHeat>`)
        - NIST interpolation for water (:ref:`waterLatentHeat.H <waterLatentHeat>`)
        - Use value for *thermophysicalProperties* dictionary (:ref:`fromThermophysicalPropertiesLatentHeat.H <fromThermophysicalPropertiesLatentHeat>`)
    - Saturation temperature/pressure (:ref:`saturationModel.H <saturationModel>`)
        - Browning Potter for sodium (:ref:`BrowningPotterSaturation.H <BrowningPotterSaturation>`)
        - NIST interpolation for water (:ref:`waterSaturation.H <waterSaturation>`)
        - TRACE model interpolation for water - YET TO BE VERIFIED (:ref:`waterTRACESaturation.H <waterTRACESaturation>`)
        - Constant temperature (see :ref:`constantTemperatureSaturation.H <constantTemperatureSaturation>`)


Regime maps
-----------

In GeN-Foam, it is possible to employ 1- and 2-dimensional regime maps to use
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

.. note::

    Anisotropic pressure drops can be set by by setting 3 different correlations
    for the 3 different local axis. An example of usage can be found in
    :ref:`FSDragFactor.H <FSDragFactor>`. In addition, it is possible to use
    anisotropic hydraulic diameter. The anisotropy of the hydraulic diameter can
    be set using the keyword *localDhAnisotropy* and assigned to it a vector of
    3 scaling factors (one for each local direction).

.. note::

    The thermal-hydraulic class can make use of a local coordinate system, which
    can be used by setting the keywords *localX* and *localY*  in the
    sub-dictionary *dragModels.(nameOfPhase).structure.(nameOfCellZones)* of the
    dictionary *constant/fluidRegion/phaseProperties*. A local coordinate system
    can be used for instance when one knows the pressure drop correlation in a
    direction that is different from the x, y, and z directions of the global
    coordinate system. Besides drag models, the local coordinate system can be
    used also for defining a tortuosity (keyword *localTortuosity*, to be
    defined as a vector in the local coordinate system).



.. toctree::
    :maxdepth: 3

    dragModels