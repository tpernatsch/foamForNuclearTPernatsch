# Material Properties and Models {#materialProperties}

The `materialProperties` keyword in the *solverDict* dictionary selects the overall treatement of material properties in OFFBEAT. Currently the code supports only the option:
- `byZone` ( see materialsByZone.H ), i.e. the cellZones (an OpenFOAM term used to indicate a subset of cells grouped under a common name) are identified with separate materials. 

It is possible for different cellZones to have the same material model, e.g. two or more cellZones can have the models and correlations characteristics of UO2. This can be useful when, for example, these cellZones differ in geometry (e.g. hollow and whole pellets in the same rod) or in some material properties (e.g. different enrichment, grain size, specific correlations for conductivity etc.).

The materials/cellZones are listed in the `materials` subdictionary in the *solverDict* file: each cellZone should have exactly one entry in the subdictionary, otherwise OFFBEAT results in an error.

OFFBEAT is provided with material models which include a set of default material properties and models. The available material models are listed [here](@ref materials).


## Thermomechanical properties

The complete list of material properties available in OFFBEAT is given hereafter.

### Conductivity models

- conductivityConstant.H , for user-selected constant value
- conductivityMatproUO2.H
- conductivityMolybdenum.H
- conductivityNfirUO2.H
- conductivityRelapZircaloy.H
- conductivityTobbe1515Ti.H
- conductivityUO2IFA601.H
- conductivityUPuO2.H

### Density models

- densityConstant.H, for user-selected constant value
- constantDensityMolybdenum.H
- constantDensityUO2.H
- constantDensityUPuO2.H
- densityIaeaZircaloy.H
- densitySchumann1515Ti.H

### Emissivity models

- emissivityConstant.H , for user-selected constant value
- constantEmissivityMolybdenum.H
- constantEmissivityZircaloy.H
- emissivityRelapUO2.H

### Heat Capacity models

- heatCapacityConstant.H , for user-selected constant value
- heatCapacityBanerjee1515Ti.H
- heatCapacityIaeaZircaloy.H
- heatCapacityMatproUO2.H
- heatCapacityMatproUPuO2.H
- heatCapacityMatproZircaloy.H
- heatCapacityMolybdenum.H

### Poisson's Ratio models

- PoissonRatioConstant.H , for user-selected constant value
- constantPoissonRatioMolybdenum.H
- constantPoissonRatioUO2.H
- constantPoissonRatioUPuO2.H
- constantPoissonRatioZircaloy.H
- PoissonRatioMatproZircaloy.H
- PoissonRatioTobbe1515Ti.H

### Thermal Expansion models

- thermalExpansionConstant.H , for user-selected constant value
- thermalExpansionGehr1515Ti.H
- thermalExpansionMatproUPuO2.H
- thermalExpansionMatproZircaloy.H
- thermalExpansionMolybdenum.H
- thermalExpansionRelapUO2.H

### Young's modulus models

- YoungModulusConstant.H , for user-selected constant value
- YoungModulusMatproUO2.H
- YoungModulusMatproUPuO2.H
- YoungModulusMatproZircaloy.H
- YoungModulusMolybdenum.H
- YoungModulusTobbe1515Ti.H

## Behavioral models

The complete list of behavioral models available in OFFBEAT is given hereafter.

### Densification models

- densificationFrapcon.H

### Failure models

- UO2meltingMatpro.H
- ZircaloyOverstrainBison.H
- ZircaloyPlasticInstabilityBison.H

### Phase Transition models

- phaseTransitionZircaloyDynamic.H

### Relocation models

- relocationFrapcon.H

### Swelling models

- swellingFrapcon.H
- swellingFrCrAl.H
- swellingGrowthAIM11515Ti.H
- swellingGrowthBisonZircaloy.H
- swellingGrowthGeneralized1515Ti.H
- swellingGrowthMatproZircaloy.H



***
Return to [Setting the `solverDict`](@ref solverDict)