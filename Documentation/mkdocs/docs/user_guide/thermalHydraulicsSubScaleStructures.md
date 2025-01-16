# Sub-scale structures

Dedicated models for specific sub-scale structures are

## Power models

Power models, i.e., active media that can provide and subtract energy, including:

- Fixed (possibly time-dependent) power (see [*fixedPower.H*](../classes/thermalHydraulics/src/phaseModels/structureModels/powerModels/fixedPower/fixedPower.md) and the tutorial *1D_CHF/imposedPower*)
- Fixed (possibly time-dependent) temperature (see [*fixedTemperature.H*](../classes/thermalHydraulics/src/phaseModels/structureModels/powerModels/fixedTemperature/fixedTemperature.md) and the tutorial *1D_CHF/imposedTemperature*)
- Heated pin, typically used for electrically heated pins (see [*heatedPin.H*](../classes/thermalHydraulics/src/phaseModels/structureModels/powerModels/heatedPin/heatedPin.md) and the tutorial *2D_KNS37-L22*)
- Nuclear fuel from FMU(s) (see [*nuclearFuelFMU.H*](../classes/thermalHydraulics/src/phaseModels/structureModels/powerModels/nuclearFuelFMU/nuclearFuelFMU.md))
- Nuclear fuel pin (see [*nuclearFuelPin.H*](../classes/thermalHydraulics/src/phaseModels/structureModels/powerModels/nuclearFuelPin/nuclearFuelPin.md) and the tutorials *3D_SmallESFR* and *2D_FFTF*)
- Lumped-parameter nuclear structure (see [*lumpedNuclearStructure.H*](../classes/thermalHydraulics/src/phaseModels/structureModels/powerModels/lumpedNuclearStructure/lumpedNuclearStructure.md)) and the tutorial *1D_thermalMSR_pointKinetics*
- Steady-state model purpose made for pebble bed reactors (see [*nuclearSteadyStatePebble.H*](../classes/thermalHydraulics/src/phaseModels/structureModels/powerModels/nuclearSteadyStatePebble/nuclearSteadyStatePebble.md)) and the tutorial *3D_gFHR*

## Heat exchanger

A heat exchanger model that is used to model the heat transfer between two disconnected regions, for instance representing the primary and secondary circuit (see *heatExchanger.H* and the tutorials *1D_HX* and *2D_FFTF*)


## Pump

A pump model used to set a (possibly time-dependent) momentum source (see *pump.H* and tutorials *2D_FFTF* and *2D_MSFR*).
