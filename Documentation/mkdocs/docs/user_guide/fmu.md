# FMU coupling

GeN-Foam provides several interface points to communicate with [Functional Mock-up Units](https://fmi-standard.org/) (FMUs). FMUs are containers of software and data that are based on a widely employed communication standard called Functional Mockup Interface (FMI). The FMI is developed by an industrial consortium led by the Modelica Association.


## Installation from source

Follow the instructions for **ECI4FOAM** (https://gitlab.com/foam-for-nuclear/ECI4FOAM) and then **FMU4FOAM** (https://gitlab.com/foam-for-nuclear/FMU4FOAM).

Here is an overview of the packages to install:
```bash
# For ECI4FOAM
sudo apt-get install libsodium-dev pybind11-dev nlohmann-json3-dev libzmq3-dev python3-dev


# For FMU4FOAM
sudo apt-get install libpugixml-dev libarchive-dev

# Python packages
python3 -m pip install omsimulator fmpy oftest pyzmq pythonfmu
```

To include the FMI features in GeN-Foam, you have to export the `LIB_ECI4FOAM` environment variable as absolute paths:
```bash
# In your ~/.bashrc, must end with "ECI4FOAM"
export LIB_ECI4FOAM="/home/.../path/to/ECI4FOAM"
export LIB_FMU4FOAM="/home/.../path/to/FMU4FOAM"

# In your ~/.bashrc, export variables for Python.h, change Python version if needed
export CPATH=/usr/include/python3.10:$CPATH
export LD_LIBRARY_PATH=/usr/lib:$LD_LIBRARY_PATH
```

Then build ECI4FOAM and FMU4FOAM. One can check that ECI4FOAM is correclty installed by running at the root of ECI4FOAM:
```bash
cd FMU4FOAM/ECI4FOAM
./Allwmake -j<n>

# Test ECI4FOAM installation
pytest

cd FMU4FOAM
./Allwmake -j<n>


```

Then build/rebuild GeN-Foam with the `--fmi` flag:
```bash
./Allwmake --fmi
# or
./Allwmake --fmi -j<N>
```


## Run a case

We recommend to use the `OMSimulator` or `FMPy` Python packages to run OpenFOAM cases with an FMU. They can be installed using the `pip` command:

```bash
pip3 install omsimulator
# and/or
pip3 install fmpy
```

Then try the [`powerTemperatureMomentumControl`](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/fmuCases/powerTemperatureMomentumControl?ref_type=heads) case to check the installation.

It is also possible to test the FMI installation using the Foam-For-Nuclear Python API on this [case](https://gitlab.com/foam-for-nuclear/foam-for-nuclear-repository/-/tree/python-api/examples/fmuCases/2D_PKCoupleFMI?ref_type=heads). This case demonstrate the use of the API to setup a case with an FMU using FMPy.


## Troubleshooting

- If your simulation crashes because of an unexpected error from externalComm json. Make sure that you have provided all the FMI ports in both codes with the correct spelling.
- If the installation of python packages is not possible, downgrade `setuptools` to `pip3 install --user --force-reinstall setuptools==65.5.1`
- Be sure to change the Python version in [`FMU4FOAM/ECI4FOAM/src/embeddingPython/Make/options`](FMU4FOAM/ECI4FOAM/src/embeddingPython/Make/options), the current version is Python3.10 but can be changed to your version.
- It is not an issue if `catch2` is not installed, this package is only present for testing and does not impact the FMI capabilities
- If `python3-dev` is not available in `apt` install, install a specific python version such as `ap install python3.10-dev`
- If `libzmq` and `cppzmq` are not installed on the machine and are not available using system installation such as `apt`, detailed explanation are provided in https://gitlab.com/foam-for-nuclear/ECI4FOAM to install these packages.
- We do not recommend to install the Python package pyfmi due to the many difficulties to install the latest version. FMPy can be used as a replacement.


## Usage using the Foam-For-Nuclear Python API

The Foam-For-Nuclear Python API can be found in this repository: https://gitlab.com/foam-for-nuclear/foam-for-nuclear-repository/-/tree/python-api?ref_type=heads. Once clone, the API can be installed using regular `pip install` in the root folder: `pip3 install .`.

This example is an extract from [2D_PKCoupleFMI](https://gitlab.com/foam-for-nuclear/foam-for-nuclear-repository/-/tree/python-api/examples/fmuCases/2D_PKCoupleFMI?ref_type=heads).

```python
# In main.py

# Coupling interface
externalCouplingDict = ffn.ExternalCouplingDict()

# Compute the total power in zone0 and output as an FMI port named gfPower
fieldIntegralToFMU = ffn.FieldIntegralToFMU(
    name="fieldIntegralToFMU",
    nameFMU="gfPower",
    fieldName="powerDensity",
    region="neutroRegion",
    cellZone="zone0",
)
externalCouplingDict.append(fieldIntegralToFMU)

# FMU simulator
FMUSimulator = ffn.FMUSimulator(
    name="FMUSimulator",
    pyClassName="ExternalReactivityController",
    pyFileName="ExternalReactivityController",
    mapping=[
        # Dir,   GeN-Foam      Modelica
        ("to",   "gfExtReact", "externalReactivity"),
        ("from", "gfPower",    "power")
    ]
)
FMUSimulator.plot_coupling_graph()

# Main model containing the case
model = ffn.Model(
    solvers=solvers,
    coupling=coupling,
    timeFolders=[timeFolder0],
    externalCouplingDict=externalCouplingDict
)
model.add_function_object(FMUSimulator)
model.export_to_openfoam()
```

```python
# In ExternalReactivityController.py

from foamForNuclear.fmi import FMPyContainer

class ExternalReactivityController(FMPyContainer):
    def __init__(
            self,
            endTime,
            jsonFile
        ):
        super().__init__(
            endTime,
            jsonFile,
            fmuName='ExternalReactivityController.fmu',
            outputFilename='ExternalReactivityController.csv'
        )
```



## Features

In this section, a list of FMI inputs/outputs in GeN-Foam is provided with a link to each class.

### Inputs from FMUs

| Feature | Location |
|:--------|:---------|
| Momentum source | [src/classes/thermalHydraulics/src/phaseModels/structureModels/pump](../classes/thermalHydraulics/src/phaseModels/structureModels/pump) |
| Fix-temperature structure | [src/classes/thermalHydraulics/src/phaseModels/structureModels/powerModels/fixedTemperature](../classes/thermalHydraulics/src/phaseModels/structureModels/powerModels/fixedTemperature/fixedTemperature.md) |
| Fix-power structure | [src/classes/thermalHydraulics/src/phaseModels/structureModels/powerModels/fixedPower](../classes/thermalHydraulics/src/phaseModels/structureModels/powerModels/fixedPower/fixedPower.md) |
| Additional heat transfer coefficient in series | [src/classes/thermalHydraulics/src/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/NusseltWallAndHfromFMU](../classes/thermalHydraulics/src/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/NusseltWallAndHfromFMU/NusseltWallAndHfromFMUFSHeatTransferCoefficient.md) |
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



---

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