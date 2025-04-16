# Porous-medium properties

The various parameters to be used in a porous-medium simulation can be set using the *phaseProperties* dictionary.


## The *phaseProperties* dictionary

The *phaseProperties* dictionary can be found in *constant/fluidRegion/*. It is a large dictionary that can be used to: choose the sub-solver to be used (one-phase, legacy one-phase or two-phase); set various properties of the phases (besides basic thermo-physical properties defined in the *thermophysicalProperties* dictionary); set the properties of the sub-scale structures (fuel pins, heat exchangers, etc) in the porous zones, including the possibility to assign a *powerModel* for power production (e.g., nuclear fuel, or constant power) and the *passiveProperties* of another sub-structure that interacts thermally with the fluid (for instance the wrappers in sodium fast reactors). In addition, models are available to model pumps and heat exchangers. The name of the porous zones must coincide with that of the cellZones of the fluidRegion mesh.

One can find detailed, commented examples in the tutorials
[3D_SmallESFR](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/reactorCases/3D_SmallESFR_NewSolverVerification/newSolver/constant/fluidRegion/phaseProperties) (single phase) and
[1D_boiling](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/featureCases/1D_boiling/constant/fluidRegion/phaseProperties) (two phases). In addition, an example of how to use a two-dimensional flow-regime map can be found in [1D_PSBT_SC](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/featureCases/1D_PSBT_SC/Phase_Ex1_12223/constant/fluidRegion/phaseProperties).


### Drag models

Currently, available models to describe pressure drops induced by the sub-scale structure or by a second phase include:

- Fluid-fluid drag models (see *FFDragCoefficientModel.H*)
	- Autruffe (see *AutruffeFFDragCoefficient.H*)
	- Bestion (see *BestionFFDragCoefficient.H*)
	- Bestion as in TRACE (see *BestionTRACEFFDragCoefficient.H*)
	- No Kazimi (see *NoKazimiFFDragCoefficient.H*)
	- Schiller Naumann (see *SchillerNaumannFFDragCoefficient.H*)
	- Wallis (see *WallisFFDragCoefficient.H*)
- Fluid-structure drag models (see *FSDragCoefficientModel.H*)
	- Baxi Dalle Donne (see *BaxiDalleDonneFSDragCoefficient.H*)
	- Churchill (see *ChurchillFSDragCoefficient.H*)
	- Engel as in TRACE (see *EngelFSDragCoefficient.H*)
	- Modified Engel (see *modifiedEngelFSDragCoefficient.H*)
	- No Kazimi (see *NoKazimiFSDragCoefficient.H*)
	- Rehme (see *RehmeFSDragCoefficient.H*)
	- Drag coefficient as Reynolds power (see *ReynoldsPowerFSDragCoefficient.H*) $A \times Re^B+C$
	- Drag coefficient as Colebrook correlation (see *ColebrookFSDragCoefficient.H*) $(A \log_{10}(Re)+B)^C$
- Two-phase drag multipliers	(see *twoPhaseDragMultiplierModel.H*)
	- Chen Kalish (see *ChenKalishTwoPhaseDragMultiplier.H*)
	- Constant (see *constantTwoPhaseDragMultiplier.H*)
	- Kaiser 74 (see *Kaiser74TwoPhaseDragMultiplier.H*)
	- Kaiser 88 (see *Kaiser88TwoPhaseDragMultiplier.H*)
	- Kottowski Savatteri (see *KottowskiSavatteriTwoPhaseDragMultiplier.H*)
	- Lockhart Martinelli (see *LockhartMartinelli.H*)
	- Lottes Flinn (see *LottesFlinnTwoPhaseDragMultiplier.H*)
	- Lottes Flinn Nguyen(see *LottesFlinnNguyenTwoPhaseDragMultiplier.H*)


### Heat transfer models

Currently, available models to describe the energy exchange with a sub-scale structure or with a second phase include:

- Fluid-fluid heat-transfer models (see *FFHeatTransferCoefficientModel.H*)
	- No Kazimi (see *NoKazimiFFHeatTransferCoefficient.H*)
	- Nusselt number correlation as Reynolds and Prandtl powers (see *NusseltFFHeatTransferCoefficient.H*) $Nu = A + B \times Re^C Pr^D$
- Fluid-structure heat-transfer models (see FSHeatTransferCoefficientModel.H)
	- Nusselt number correlation as Reynolds and Prandtl powers and surface to bulk temperature ratio (see *NusseltFSHeatTransferCoefficient.H*) $Nu = A + B \times Re^C Pr^D \left( \frac{T_w}{T_b} \right)^E$
	- Nusselt number correlation as Reynolds and Prandtl powers, plus an additional heat transfer coefficient to take into account the resistance of a wall (see *NusseltAndWallFSHeatTransferCoefficient.H*) $H = \frac{Nu \times \kappa}{D_h} + H_{wall} \quad \text{with} \quad Nu = A + B \times Re^C Pr^D$
	-  Shah (see *ShahFSHeatTransferCoefficient.H*)
	-  Gorenflo (see *GorenfloFSHeatTransferCoefficient.H*)
	-  multiRegimeBoilingTRACE - multi-regime heat transfer coefficient that replicates what TRACE does below CHF (see *multiRegimeBoilingTRACEFSHeatTransferCoefficient.H*) ADD REF GAUTHIER
	-  multiRegimeBoilingTRACE - multi-regime heat transfer coefficient that replicates what TRACE does, including CHF and post-CHF. Not verified! It requires specifying the *multiRegimeBoilingTRACECHF* model for the water.structure heat transfer, and the *multiRegimeBoilingVapourTRACE* model for the vapour.structure heat transfer (see *multiRegimeBoilingTRACECHFFSHeatTransferCoefficient.H* and *multiRegimeBoilingVapourTRACEFSHeatTransferCoefficient.H*). Please notice that a lookup table for CHF is still missing.
	-  multiRegimeBoiling - multi-regime heat transfer coefficient that replicates the same logic as TRACE, but with more flexibility for user-selectable sub-models (see *multiRegimeBoiling.H*). Can be used below CHF.
	-  Sub-models employed by the multi-regime models:
		- Critical heat flux related models (see *CHFModel.H*):
			- Critical heat flux models
				- Constant, user-selectable value (see *constantCHF.H*)
				- Lookup table, not yet implemented (empty class at *lookUpTableCHF.H*)
			- Leidenfrost models (see *TLFModel.H*)
				- Groeneveld Stewart (see *GroeneveldStewartTLF.H*)
		- Flow Enhancement Factor Models (see *flowEnhancementFactorModel.H*)
			- Chen (see *ChenFlowEnhancementFactor.H*)
			- COBRA-TF (see *COBRA-TFFlowEnhancementFactor.H*)
			- Rezkallah Sims (see *RezkallahSimsFlowEnhancementFactor.H*)
		- Post-CHF models
			- Cachard (for liquid) (see *CachardLiquidFSHeatTransferCoefficient.H*)
			- Cachard (for vapour) (see *CachardVapourFSHeatTransferCoefficient.H*)
		- Sub-Cooled Boiling Fraction Models (see *subCooledBoilingFractionModel.H*)
			- Constant (see *constantSubCooledBoilingFraction.H*)
			- Saha Zuber (see *SahaZuberSubCooledBoilingFraction.H*)
		- Superposition Nucleate Boiling (see *superpositionNucleateBoilingFSHeatTransferCoefficient.H*)
		- Suppression factor models (see *suppressionFactorModel.H*)
			- Chen (see *ChenSuppressionFactor.H*)
			- COBRA-TF (see *COBRA-TFSuppressionFactor.H*)
		- Temperature of the onset of nucleate boiling (see *TONBModel.H*)
			- Basu (see *BasuTONB.H*)


### Special models

Dedicated models for specific sub-scale structures are described in more details in [this page](thermalHydraulicsSubScaleStructures.md).


### Models for two-phase flows

Currently, available models for two-phase flow simulations include:

- Contact partition models (see *contactPartitionModel.H*)
	- Linear (see *linearContactPartition.H*)
	- Complementary (see *complementaryContactPartition.H*)
- Dispersions models (see *dispersionModel.H*)
	- Constant (see *constantDispersion.H*)
- Fluid diameter models (see *fluidDiameterModel.H*)
	- Iso-molar bubble (see *isomolarBubbleFluidDiameter.H*)
	- Iso-thermal bubble  (see *isothermalBubbleFluidDiameter.H*)
	- Pipe film (see *pipeFilmFluidDiameter.H*)
	- Wallis film (see *WallisFilmFluidDiameter.H*)
- Virtual mass models (see *virtualMass.H*)
	- Virtual mass coefficient (see *virtualMassCoefficientModel.H*)
- Interfacial area models (see *interfacialAreaModel.H*)
	- Annular (see *annularInterfacialArea.H*)
	- No Kazimi (see *NoKazimiInterfacialArea.H*)
	- Schor (see *SchorInterfacialArea.H*)
	- Spherical (see *sphericalInterfacialArea.H*)
- Phase change models
	- Forced constant (see *forcedConstantPhaseChange.H*)
	- Heat driven (see *heatDrivenPhaseChange.H*)
	- Latent heat (see *latentHeatModel.H*)
		- Fink Leibowitz for sodium (see *FinkLeibowitzLatentHeat.H*)
		- NIST interpolation for water (see *waterLatentHeat.H*)
		- Use value for *thermophysicalProperties* dictionary (see *fromThermophysicalPropertiesLatentHeat.H*)
	- Saturation temperature/pressure (see *saturationModel.H*)
		- Browning Potter for sodium (see *BrowningPotterSaturation.H*)
		- NIST interpolation for water (see *waterSaturation.H*)
		- TRACE model interpolation for water - YET TO BE VERIFIED (see *waterTRACESaturation.H*)
		- Constant temperature (see *constantTemperatureSaturation.H*)


### Regime maps

In GeN-Foam, it is possible to employ 1- and 2-dimensional regime maps to use different models for different flow conditions. This can be used for instance in one-phase simulation to provide different correlations for turbulent and laminar flow (see [3D_SmallESFR](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/reactorCases/3D_SmallESFR_NewSolverVerification/newSolver/constant/fluidRegion/phaseProperties) for a commented example), or in 2-phase flow simulations to provide full regime maps (see [1D_boiling](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/featureCases/1D_boiling/constant/fluidRegion/phaseProperties) for a commented example of a 1-dimensional map, and [1D_PSBT_SC](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/featureCases/1D_PSBT_SC/Phase_Ex1_12223/constant/fluidRegion/phaseProperties) for a commented example of a 2-dimensional map). Multiple maps can be used in the same simulation. In addition, regime-maps models can be mixed with multi-regime models, as in [1D_PSBT_SC](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/featureCases/1D_PSBT_SC/Phase_Ex1_12223/constant/fluidRegion/phaseProperties), where a *preCHFTraceRegimeMap* is employed to assign models for phase dispersion, interfacial area and bubble diameter, while a single multi-regime model is employed to describe heat transfer between liquid and structure throughout the various regimes.

N.B.: Anisotropic pressure drops can be set by by setting 3 different correlations for the 3 different local axis. An example of usage can be found in FSDragFactor.H . In addition, it is possible to use anisotropic hydraulic diameter. The anisotropy of the hydraulic diameter can be set using the keyword *localDhAnisotropy* and assigned to it a vector of 3 scaling factors (one for each local direction).

N.B.2: The thermal-hydraulic class can make use of a local coordinate system, which can be used by setting the keywords *localX* and *localY*  in the sub-dictionary *dragModels.(nameOfPhase).structure.(nameOfCellZones)* of the dictionary *constant/fluidRegion/phaseProperties*. A local coordinate system can be used for instance when one knows the pressure drop correlation in a direction that is different from the x, y, and z directions of the global coordinate system. Besides drag models, the local coordinate system can be used also for defining a tortuosity (keyword *localTortuosity*, to be defined as a vector in the local coordinate system).
