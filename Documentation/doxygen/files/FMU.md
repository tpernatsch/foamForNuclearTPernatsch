# FMU coupling {#FMU}

GeN-Foam provides several interface points to communicate with [Functional Mock-up Units](https://fmi-standard.org/) (FMUs). FMUs are containers of software and data that are based on a widely employed communication standard called Functional Mockup Interface (FMI). The FMI is developed by an industrial consortium led by the Modelica Association.



## Compiling

To use the FMI coupling interface in GeN-Foam, the user has to install the [FMU4FOAM](https://github.com/DLR-RY/FMU4FOAM) project developed by the DLR using the following commands:

First of all, you have to install an old version of **conan**. Nothing works with 2.x
```
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

To include into GeN-Foam you have to uncomment the commented lines in GeN-Foam/Make/options. Then the GeN-Foam project can be built as usual.


## Features

In this section, a list of FMI inputs/outputs in GeN-Foam is provided with a link to each class.

### Inputs from FMUs

| Feature | Location |
|:--------|:---------|
| Momentum source | [GeN-Foam/classes/thermalHydraulics/src/phaseModels/structureModels/pump](GeN-Foam/classes/thermalHydraulics/src/phaseModels/structureModels/pump) |
| Fix-temperature structure | [GeN-Foam/classes/thermalHydraulics/src/phaseModels/structureModels/powerModels/fixedTemperatureFMU](GeN-Foam/classes/thermalHydraulics/src/phaseModels/structureModels/powerModels/fixedTemperatureFMU) |
| External reactivity in the point-kinetics solver | [GeN-Foam/classes/neutronics/pointKinetics](GeN-Foam/classes/neutronics/pointKinetics) |
| External neutron source modulation in the point-kinetics solver | [GeN-Foam/classes/neutronics/pointKinetics](GeN-Foam/classes/neutronics/pointKinetics) |
| Additional heat transfer coefficient in serie | [GeN-Foam/classes/thermalHydraulics/src/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/NusseltWallAndHfromFMU](GeN-Foam/classes/thermalHydraulics/src/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/NusseltWallAndHfromFMU) |


### Outputs to FMUs:

| Feature | Location |
|:--------|:---------|
| External sensor | FMU4FOAM |
| Field integral over a cellZone | [GeN-Foam/classes/thermalHydraulics/src/functionObjects/fieldIntegralToFMU](GeN-Foam/classes/thermalHydraulics/src/functionObjects/fieldIntegralToFMU) |


## Potential Issues

Make sure that `-I./*/FMU4FOAM/ECI4FOAM/src/externalComm/lnInclude` and `-lexternalComm` have been uncommented in:
- [GeN-Foam/classes/neutronics/Make/options](GeN-Foam/classes/neutronics/Make/options)
- [GeN-Foam/classes/thermalHydraulics/src/Make/options](GeN-Foam/classes/thermalHydraulics/src/Make/options)
- [GeN-Foam/Make/options](GeN-Foam/Make/options)



