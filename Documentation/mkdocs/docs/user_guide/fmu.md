# FMU coupling

GeN-Foam provides several interface points to communicate with [Functional Mock-up Units](https://fmi-standard.org/) (FMUs). FMUs are containers of software and data that are based on a widely employed communication standard called Functional Mockup Interface (FMI). The FMI is developed by an industrial consortium led by the Modelica Association.


## Compiling

Follow the instructions for
- ECI4FOAM: https://gitlab.com/foam-for-nuclear/ECI4FOAM
- FMU4FOAM: https://gitlab.com/foam-for-nuclear/FMU4FOAM

To include into GeN-Foam you have to export the `LIB_ECI4FOAM` environment variable such as:
```bash
# In your .bashrc, must end with "ECI4FOAM"
export LIB_ECI4FOAM="/home/.../path/to/ECI4FOAM"
```

Then execute:
```bash
./Allwmake --fmi
# or
./Allwmake --fmi -j<N>
```


### Alternative installation from the original FMU4FOAM repository

To use the FMI coupling interface in GeN-Foam, the user has to install the [FMU4FOAM](https://github.com/DLR-RY/FMU4FOAM) project developed by the DLR using the following commands:

First of all, you have to install an old version of **conan**. Nothing works with 2.x
```bash
pip install oftest conan==1.58.0
```

Now clone the repository into the src folder:
```bash
cd src
git clone https://github.com/DLR-RY/FMU4FOAM.git
```

Now run:
```bash
cd FMU4FOAM
./build-ECI4FOAM.sh
```

Then run the `Allwmake` in `ECI4FOAM` and in the `FMU4FOAM` root folder:
```bash
cd ECI4FOAM
./Allwmake
# Get to the root folder of FMU4FOAM
cd ..
./Allwmake
```

Finally, install the Python packages for FMUs:
```bash
pip install fmu4foam OMSimulator pythonfmu
```

To test:
```bash
cd examples/heatedRoom
./Allrun
```

To include into GeN-Foam you have to export the `LIB_ECI4FOAM` environment variable such as:
```bash
# In your .bashrc, must end with "ECI4FOAM"
export LIB_ECI4FOAM="/home/.../path/to/ECI4FOAM"
```
Then the GeN-Foam project can be built as usual.


## Features

In this section, a list of FMI inputs/outputs in GeN-Foam is provided with a link to each class.

### Inputs from FMUs

| Feature | Location |
|:--------|:---------|
| Momentum source | [src/classes/thermalHydraulics/src/phaseModels/structureModels/pump](../classes/thermalHydraulics/src/phaseModels/structureModels/pump) |
| Fix-temperature structure | [src/classes/thermalHydraulics/src/phaseModels/structureModels/powerModels/fixedTemperature](../classes/thermalHydraulics/src/phaseModels/structureModels/powerModels/fixedTemperature/fixedTemperature.md) |
| Fix-power structure | [src/classes/thermalHydraulics/src/phaseModels/structureModels/powerModels/fixedPower](../classes/thermalHydraulics/src/phaseModels/structureModels/powerModels/fixedPower/fixedPower.md) |
| Additional heat transfer coefficient in serie | [src/classes/thermalHydraulics/src/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/NusseltWallAndHfromFMU](../classes/thermalHydraulics/src/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/NusseltWallAndHfromFMU/NusseltWallAndHfromFMUFSHeatTransferCoefficient.md) |
| External reactivity in the point-kinetics solver | [src/classes/neutronics/pointKinetics](../classes/neutronics/pointKinetics/pointKineticNeutronics.md) |
| External neutron source modulation in the point-kinetics solver | [src/classes/neutronics/pointKinetics](../classes/neutronics/pointKinetics/pointKineticNeutronics.md) |
| Boron reactivity in the point-kinetics solver | [src/classes/neutronics/pointKinetics](../classes/neutronics/pointKinetics/pointKineticNeutronics.md) |
| Decay power in the point-kinetics solver | [src/classes/neutronics/pointKinetics](../classes/neutronics/pointKinetics/pointKineticNeutronics.md) |
| Time profile object | [src/classes/common/timeProfile](../classes/common/timeProfile) |


### Outputs to FMUs

| Feature | Location |
|:--------|:---------|
| External sensor | FMU4FOAM |
| Field integral over a cellZone | [src/classes/thermalHydraulics/src/functionObjects/fieldIntegralToFMU](../classes/thermalHydraulics/src/functionObjects/fieldIntegralToFMU/fieldIntegralToFMU.md) |


### Inputs and outputs from/to an FMU

| Feature | Location |
|:--------|:---------|
| Nuclear fuel structure based on an FMU | [src/classes/thermalHydraulics/src/phaseModels/structureModels/powerModels/nuclearFuelFMU](../classes/thermalHydraulics/src/phaseModels/structureModels/powerModels/nuclearFuelFMU/nuclearFuelFMU.md) |


## Tips

If your simulation crashes because of an unexpected error from externalComm json. Make sure that you have provided all the FMI ports in both codes with the correct spelling.
