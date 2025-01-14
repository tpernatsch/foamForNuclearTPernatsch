# Thermal-hydraulics

## Introduction

Both single- and two-phase simulations can be performed using GeN-Foam. All sub-solvers were developed for a coarse-mesh porous-medium treatment of complex structures such as core and heat exchanger, and for a standard RANS treatment of clear-fluid regions. The sub-solvers automatically switch from a porous-medium (coarse-mesh) treatment to a standard CFD (fine-mesh) treatment when the volume fraction of the sub-scale structures is set to zero. This allows for an implicit coupling of porous-medium (sub-channel-like in 2D and 3D, or system-code-like) treatment of complex structures (e.g., core and heat exchangers) with a standard CFD treatment of clear-fluid regions (e.g., plena and pools).

A coarse-mesh porous-medium treatment of the core implies that the core is modeled without resolving the sub-scale structure (e.g., the fuel rods or the heat exchanger tubes). As a matter of fact, in principle and for consistency, the finest radial mesh chosen by a user should not be finer than one cell per pin cell. A porous-medium formulation derives from a volume averaging of the Navier-Stokes equations. The volume averaging results in source terms that describe the interaction (drag and heat transfer) of the fluid with the sub-scale structure. In GeN-Foam, these source terms are modeled using user-selectable correlations for drag (e.g., correlations for the Darcy friction factor) and heat transfer (e.g., correlations for the Nusselt number). In this sense, a porous-medium model can be associated with a 3-D version of a system code.

With regards to the modeling of the sub-scale structures, GeN-Foam allows modeling simultaneously in the same region both a "power model" and a "passive structure". Power models are used to model for instance the nuclear fuel (based on a 1-D approximation), electrically heated rods, or a fixed temperature body (which can be used to approximate a heat exchanger). Passive structures are structures that passively heat up or cool down based on their heat capacity, volumetric area, and heat transfer with the coolant. This can be used to model structures like assembly wrappers or reflectors.

All thermal-hydraulics functionalities are handled by the class *thermalHydraulicsModel.H*, the derived classes for the various sub-solvers (see below), and a thermal-hydraulic library that can be found under */GeN-Foam/classes/thermalHydraulics/src*.


### The porous-medium approach in GeN-Foam

GeN-Foam was born for safety analyses and, to reduce computational footprint, its base approach is to model for instance the core as a porous medium. In a porous-medium approach, the fuel and other structures (e.g., the assembly wrappers) are modeled using sub-scale models. This means that, in each cell, we have the fluid and one or two lumped models for the sub-scale structures. The simplest structures in GeN-Foam are the passive structures. These passive structures are simply modeled as a heat capacity and they can be used for instance to model assembly wrappers or reflector structures. In essence, the fluid will interact cell-by-cell with these passive structures: they will take energy from the fluid if the temperature of the fluid is higher than that of the surface of the structure, and vice versa. A more complex example of sub-scale structure is given by the powerModels. An example of a power model is the nuclearFuelPin, which can be used to model a standard pin-type fuel. This power model is capable of getting the power density from neutronics, solving a 1-D model for heat transfer in the fuel, and giving back to the fluid the temperature at the surface of the cladding. The fluid will then be capable of calculating the heat transfer with the fuel based on the cladding surface temperature and the Nusselt number. Each cellZone can host one passive structure and one powerModel.


## Sub-solvers

Thermal-hydraulics calculations are performed by classes derived from *thermalHydraulicsModel.H* that contain specific sub-solvers:
- *onePhase* for single-phase calculations, using the formulation proposed in Refs. [@Radman2019ADesign] [@RADMAN2021111178] [@RADMAN2021111422]  (see *onePhase.H*)
- *twoPhase* for adjoint diffusion calculations, using the formulation proposed in Refs. [@Radman2019ADesign] [@RADMAN2021111178] [@RADMAN2021111422] (see *twoPhase.H*)

For the user, the derived classes translate into runtime selectable models. The specific sub-solver to be used in a simulation can be selected in the solvers dictionary like explained in the coupling section.


### OpenFOAM-based sub-solvers

The new GeN-Foam structure allows to include already developed OpenFOAM solver. This requires to transpose the solver application into a GeN-Foam solver class. *compressibleInterFoam* shows how to translate standard OpenFOAM solvers into this new format, taking into account all the dependencies.

We list below the imported OpenFOAM-based standard solvers available in GeN-Foam.
* *compressibleInterFoam* for two compressible, non-isothermal immiscible fluids using a VOF (volume of fluid) phase-fraction based interface capturing approach. This standard solver has been transposed from OpenFOAM ([link](https://www.openfoam.com/documentation/guides/latest/man/compressibleInterFoam.html)) into the new GeN-Foam solver structure (see *src/classes/openFoamImportedSolvers/compressibleInterFoam*).


## Porous-medium properties

The various parameters to be used in a porous-medium simulation can be set using the *phaseProperties* dictionary.


### The *phaseProperties* dictionary

The *phaseProperties* dictionary can be found in *constant/fluidRegion/*. It is a large dictionary that can be used to: choose the sub-solver to be used (one-phase, legacy one-phase or two-phase); set various properties of the phases (besides basic thermo-physical properties defined in the *thermophysicalProperties* dictionary); set the properties of the sub-scale structures (fuel pins, heat exchangers, etc) in the porous zones, including the possibility to assign a *powerModel* for power production (e.g., nuclear fuel, or constant power) and the *passiveProperties* of another sub-structure that interacts thermally with the fluid (for instance the wrappers in sodium fast reactors). In addition, models are available to model pumps and heat exchangers. The name of the porous zones must coincide with that of the cellZones of the fluidRegion mesh.

One can find detailed, commented examples in the tutorials
[3D_SmallESFR](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/reactorCases/3D_SmallESFR_NewSolverVerification/newSolver/constant/fluidRegion/phaseProperties) (single phase) and
[1D_boiling](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/featureCases/1D_boiling/constant/fluidRegion/phaseProperties) (two phases). In addition, an example of how to use a two-dimensional flow-regime map can be found in [1D_PSBT_SC](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/featureCases/1D_PSBT_SC/Phase_Ex1_12223/constant/fluidRegion/phaseProperties).

**Drag models**. Currently, available models to describe pressure drops induced by the sub-scale structure or by a second phase include:

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

**Heat transfer models**. Currently, available models to describe the energy exchange with a sub-scale structure or with a second phase include:

- Fluid-fluid heat-tranfer models (see *FFHeatTransferCoefficientModel.H*)
	- No Kazimi (see *NoKazimiFFHeatTransferCoefficient.H*)
	- Nusselt number correlation as Reynolds and Prandtl powers (see *NusseltFFHeatTransferCoefficient.H*) $Nu = A + B \times Re^C Pr^D$
- Fluid-structure heat-tranfer models (see FSHeatTransferCoefficientModel.H)
	- Nusselt number correlation as Reynolds and Prandtl powers and surface to bulk temperature ratio (see *NusseltFSHeatTransferCoefficient.H*) $Nu = A + B \times Re^C Pr^D \left( \frac{T_w}{T_b} \right)^E$
	- Nusselt number correlation as Reynolds and Prandtl powers, plus an additional heat transfer coefficient to take into account the resistance of a wall (see *NusseltAndWallFSHeatTransferCoefficient.H*) $H = \frac{Nu \times \kappa}{D_h} + H_{wall} \quad \text{with} \quad Nu = A + B \times Re^C Pr^D$
	-  Shah (see *ShahFSHeatTransferCoefficient.H*)
	-  Gorenflo (see *GorenfloFSHeatTransferCoefficient.H*)
	-  multiRegimeBoilingTRACE - multi-regime heat transfer coefficient that replicates what TRACE does below CHF (see *multiRegimeBoilingTRACEFSHeatTransferCoefficient.H*) ADD REF GAUTHIER
	-  multiRegimeBoilingTRACE - multi-regime heat transfer coefficient that replicates what TRACE does, including CHF and post-CHF. Not verified! It requires specifying the *multiRegimeBoilingTRACECHF* model for the water.structure heat transfer, and the *multiRegimeBoilingVapourTRACE* model for the vapour.structure heat tranfer (see *multiRegimeBoilingTRACECHFFSHeatTransferCoefficient.H* and *multiRegimeBoilingVapourTRACEFSHeatTransferCoefficient.H*). Please notice that a lookup table for CHF is still missing.
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

**Special models**. Dedicated models for specific sub-scale structures include:

- Power models, i.e., active media that can provide and subtract energy, including:
	- Fixed (possibly time-dependent) power (see *fixedPower.H* and the tutorial *1D_CHF/imposedPower*)
	- Fixed (possibly time-dependent) temperature (see *fixedTemperature.H* and the tutorial *1D_CHF/imposedTemperature*)
	- Heated pin, typically used for electrically heated pins (see *heatedPin.H* and the tutorial *2D_KNS37-L22*)
	- Nuclear fuel from FMU(s) (see *nuclearFuelFMU.H*)
	- Nuclear fuel pin (see *nuclearFuelPin.H* and the tutorials *3D_SmallESFR* and *2D_FFTF*)
	- Lumped-parameter nuclear structure (see *lumpedNuclearStructure.H*) and the tutorial *1D_thermalMSR_pointKinetics*
	- Steady-state model purpose made for pebble bed reactors (see *nuclearSteadyStatePebble.H*) and the tutorial *3D_gFHR*
- A heat exchanger model that is used to model the heat transfer between two disconnected regions, for instance representing the primary and secondary circuit (see *heatExchanger.H* and the tutorials *1D_HX* and *2D_FFTF*)
- A pump model used to set a (possibly time-dependent) momentum source (see *pump.H* and tutorials *2D_FFTF* and *2D_MSFR*).

**Models for two-phase flows**. Currently, available models for two-phase flow simulations include:

- Contact partition models (see *contactPartitionModel.H*)
	- Linear (see *linearContactPartition.H*)
	- Complementary (see *complementaryContactPartition.H*)
- Disperions models (see *dispersionModel.H*)
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
	- No Kazimi  (see *NoKazimiInterfacialArea.H*)
	- Schor (see *SchorInterfacialArea.H*)
	- Spherical (see *sphericalInterfacialArea.H*)
- Phase change models
	- Forced constant (see *forcedConstantPhaseChange.H*)
	- Heat driven  (see *heatDrivenPhaseChange.H*)
	- Latent heat (see *latentHeatModel.H*)
		- Fink Leibowitz for sodium (see *FinkLeibowitzLatentHeat.H*)
		- NIST interpolation for water (see *waterLatentHeat.H*)
		- Use value for thermophysicalProerties dictionary (see *fromThermophysicalPropertiesLatentHeat.H*)
	- Saturation temperature/pressure (see *saturationModel.H*)
		- Browning Potter for sodium (see *BrowningPotterSaturation.H*)
		- NIST interpolation for water (see *waterSaturation.H*)
		- TRACE model interpolation for water - YET TO BE VEIFIED (see *waterTRACESaturation.H*)
		- Constant temperature (see *constantTemperatureSaturation.H*)

**Regime maps**. In GeN-Foam, it is possible to employ 1- and 2-dimensional regime maps to use different models for different flow conditions. This can be used for instance in one-phase simulation to provide different correlations for turbulent and laminar flow (see [3D_SmallESFR](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/reactorCases/3D_SmallESFR_NewSolverVerification/newSolver/constant/fluidRegion/phaseProperties) for a commented example), or in 2-phase flow simulations to provide full regime maps (see [1D_boiling](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/featureCases/1D_boiling/constant/fluidRegion/phaseProperties) for a commented example of a 1-dimensional map, and [1D_PSBT_SC](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/featureCases/1D_PSBT_SC/Phase_Ex1_12223/constant/fluidRegion/phaseProperties) for a commented example of a 2-dimensional map). Multiple maps can be used in the same simulation. In addition, regime-maps models can be mixed with multi-regime models, as in [1D_PSBT_SC](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/featureCases/1D_PSBT_SC/Phase_Ex1_12223/constant/fluidRegion/phaseProperties), where a *preCHFTraceRegimeMap* is employed to assign models for phase dispersion, interfacial area and bubble diameter, while a single multi-regime model is employed to describe heat transfer between liquid and structure throughout the various regimes.

N.B.: Anisotropic pressure drops can be set by by setting 3 different correlations for the 3 different local axis. An example of usage can be found in FSDragFactor.H . In addition, it is possible to use anisotropic hydraulic diameter. The anisotropy of the hydraulic diameter can be set using the keyword *localDhAnisotrpy* and assigned to it a vector of 3 scaling factors (one for each local direction).

N.B.2: The thermal-hydraulic class can make use of a local coordinate system, which can be used by setting the keywords *localX* and *localY*  in the sub-dictionary *dragModels.(nameOfPhase).structure.(nameOfCellZones)* of the dictionary *constant/fluidRegion/phaseProperties*. A local coordinate system can be used for instance when one knows the pressure drop correlation in a direction that is different from the x, y, and z directions of the global coordinate system. Besides drag models, the local coordinate system can be used also for defining a tortuosity (keyword *localTortuosity*, to be defined as a vector in the local coordinate system).


## Physical properties

### The *g* dictionary

The *g* dictionary can be found under *constant/fluidRegion/*. It is a standard OpenFOAM dictionary that allows specifying the gravitational acceleration.


### The *thermophysicalProperties* dictionary

The *thermophysicalProperties* dictionary can be found under *constant/fluidRegion/*. It is a standard OpenFOAM dictionary that allows defining the thermo-physical properties of the coolant. When performing two-phase flow analyses, two dictionaries must be employed named*thermophysicalProperties.(name of fluid)*. The name of the two fluids is defined in the *phaseProperties* dictionary.

One can find a detailed, commented example in the tutorials
[3D_SmallESFR](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/reactorCases/3D_SmallESFR_NewSolverVerification/newSolver/constant/fluidRegion/thermophysicalProperties) (one-phase),
[1D_boiling (liquid)](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/featureCases/1D_boiling/constant/fluidRegion/thermophysicalProperties.liquid) (two-phase, liquid)
[1D_boiling (vapour)](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/featureCases/1D_boiling/constant/fluidRegion/thermophysicalProperties.vapour) (two-phase, vapour)


## Turbulence properties

### The *turbulenceProperties* dictionary

The *turbulenceProperties* dictionary can be found under *constant/fluidRegion/*. It is a standard OpenFOAM dictionary that allows defining the turbulence model to be used.

When clear-fluid simulations (i.e., without porous zones) are performed, one can use the standard kEpsilon model of OpenFOAM.

When porous zones are present in the simulation, it is recommended to use *porousKEpsilon* (see *porousKEpsilon.H*). The only difference w.r.t. the standard k-epsilon model is that it forces k and epsilon to equilibrium values inside the porous zones. These equilibrium values can be set in the *porousKepsionProperties* sub-dictionary. Please notice that a porous medium simulation using the equilibrium values of k and epsilon for the sub-scale structure (viz., the values inside a fuel sub-channel) would entail the risk of an unstable solution. This occurs because the turbulent viscosity is primarily associated with the sub-scale structure, which might not be sufficient to maintain stability at the coarse mesh's length scale. To address this problem, one can define the keyword DhStruct in *constant/fluidRegion/phaseProperties/dragModels.(nameOfPhase).structure.(nameOfCellZones)*. This keyword defines the hydraulic diameter of the whole porous structure (viz., the dimension of the assembly, if using baffles to model wrappers, or of the entire core). The code uses it to make sure the turbulent viscosity results in a laminar Reynolds number (defaulted to 500).


While some approaches to model k and epsilon for two-phase flow simulations are presently included in the code. In particular, the Lahey model (see *LaheyKEpsilon.H*) and a mixture model (see *mixtureKEpsilon.H*) can be used for clear-fluids, or mixed clear-fluid and porous-medium simulations in case of strongly advective two-phase flow scenarios where turbulent mixing may be neglected. In addition, a simple extension of the *porousKEpsilon* model has been implemented that allows to correct the turbulent intensity using a term that is proportional to the fraction of the other phase (see *porousKEpsilon2PhaseCorrected.H*).

One can find a detailed, commented example of a porous one-phase simulation in the tutorial
[3D_SmallESFR](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/reactorCases/3D_SmallESFR_NewSolverVerification/newSolver/constant/fluidRegion/turbulenceProperties).


## Initial and boundary conditions

Initial and boundary conditions adopt the usual OpenFOAM logic for one- and two-phase solvers. A couple of things to be kept in mind:

- The pressure field we solve for is p_rgh (pressure minus the gravitational head)
- When performing turbulent analyses, one needs to add the fields *k*, *epsilon*, nut and *alphat*

One thing that instead specific to GeN-Foam (except for the one-phase legacy sub-solver) and that one needs to keep in mind is that U (or u.(name of fluid)) are the real velocities, not the Darcy velocities. In a porous structure, they represent the actual velocity of the fluid, and not the velocity multiplied by the fluid fraction. For instance, U will increase when transiting from a high- to a low-porosity region.

OpenFOAM provides most of the boundary conditions one may need for thermal-hydraulics models. In addition, a few boundary conditions have been included in GeN-Foam and can be found in *GeN-Foam/classes/thermalHydraulics/src/boundaryConditions*. Information on the use of each boundary condition can be found in the header files (.H).


### Setting the initial power

There are several ways to set the power in GeN-Foam. For the power generated in subscale structures:

- The thermal-hydraulics solver will normally use the *powerDensity.* fields (for instance, *powerDensity.nuclearFuelPin* for pin-based reactors) that it finds in the 0 (or *startTime*) folder.
- As an alternative, one can provide cellZone-by-cellZone values via the keyword *powerDensity* in the various power models in the *phase* properties (see for instance [1D_boiling](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/featureCases/1D_boiling/constant/fluidRegion/phaseProperties)). However, if GeN-Foam finds the corresponding field in the 0 (or *startTime*) folder, this will take priority.

For the power generated in the fluid itself:

- The thermal-hydraulics solver will normally use the *powerDensity* field that it finds in the 0 (or *startTime*) folder.
- One can override this behavior by using the *initialPowerDensity* keyword in the *phaseProperties* (see [1D_MSR_pointKinetics](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/reactorCases/1D_MSR_pointKinetics/rootCase/constant/fluidRegion/phaseProperties)) for an example. Also in this case the field in the 0 (or *startTime*) folder will take priority.

If neutronics is activated, the power density will be taken from the neutronic sub-solver. For eigenvalue calculations, the power is set in the *pTarget* keyword in the *reactorState* dictionary. For transients, the power is a result of calculations. There is one important exception to this behavior: the point kinetics solver will only rescale the power densities (see below about why plural) that it finds in *neutroRegion*, or, if it does not find it, the one(s) that it finds in *fluidRegion*. The rescaled power density will be written to both *neutroRegion* and *fluidRegion*. For point kinetics, the *pTarget* keyword in *reactorState* is not used by the solver itself. However, to correctly plot point kinetics results, pTarget must be consistent with the mentioned power densities.

NB: In two-phase simulations with liquid fuel, the powerDensity in neutronics goes to anything that is liquid in thermal-hydraulics. You are supposed to have one liquid and one gas. Otherwise, power will be counted twice.


### Power densities and secondary power densities, and liquid fuel

The spatial neutronics solvers always create a *powerDensity* and *secondaryPowerDensity* fields. By default, *secondaryPowerDensity* is set to zero and the *fuelFraction* keyword in *nuclearData* is used to translate the volume-average power density that is normally calculated by multiplying cross-sections and fluxes into the fuel-averaged power density that is needed by the thermal-hydraulic sub-solver

However, a *secondaryPowerDensity* might sometimes be needed. It might be used to provide some power to the coolant in a solid-fuel reactor and, more importantly, to provide some power to the graphite in a liquid-fuel reactor. To calculate a *secondaryPowerDensity*, GeN-Foam needs to know how much of the total power goes into the *secondaryPowerDensity*, and what is the volume fraction of the secondary power-producing structure or liquid. This can be done by using the *fractionToSecondaryPower* and *secondaryPowerVolumeFraction* keywords in each cellZone in *nuclearData* (the same place as *fuelFraction*). If these keywords are present, GeN-Foam will calculate power densities as follows:

- $\text{secondaryPowerDensity} = \frac{\text{powerDensity}}{\max(\text{secondaryPowerVolumeFraction}, \text{SMALL}) \times \text{fractionToSecondaryPower}}$
- $\text{powerDensity} = \frac{\text{powerDensity}}{\max(\text{fuelFraction}, SMALL) \times (1.0 - \text{fractionToSecondaryPower})}$

When the *liquidFuel* flag is set to false, the thermal-hydraulic sub-solver will:

- take the *powerDensity* field from neutronics and project it to its own *powerDensityNeutronics* field;
- take the *secondaryPowerDensity* field from neutronics and project it to its own *powerDensityNeutronicsToLiquid* field.

When the *liquidFuel* flag is set to true, the thermal-hydraulic sub-solver will:

- take the *powerDensity* field from neutronics and project it to its own *powerDensityNeutronicsToLiquid* field;
- take the *secondaryPowerDensity* field from neutronics and project it to its own *powerDensityNeutronics* field.

When point kinetics is used, the solver will simply rescale the *powerDensity* and *secondaryPowerDensity* it finds, and the thermal-hydraulic solver will take them depending on the *liquidFuel* flag as described above. The only exception is when *liquidFuel* is true and the *initialPowerDensity* keyword is used. In this case, *initialPowerDensity* will take priority and this is the value that GeN-Foam will rescale and print to the powerDensityToLqiuid


NB1: The power densities in the thermal-hydraulic sub-solver are ALWAYS the physical ones: for instance, when the *nuclearFuelPin* model is used for pin-based reactors, *powerDensity* refers to the power density inside the fuel matrix. For liquid fuel, the *powerDensity* is the power density in the liquid. They are not the power densities smeared over the whole volume.

NB2: If you calculate the powerDensity using Serpent, you have to divide it by the fuel fractio before feeding it to GeN-Foam (see [Neutronics](neutronics.md), in the introduction, NB2)


## Discretization and solution

Details for discretization and solution of equations are handled in a standard OpenFOAM way, i.e., through the *fvSolution* and *fvSchemes* dictionaries in *constant/fluidRegion*.
