# How to install FMU4FOAM

First of all, you have to install an old version of **conan**. Nothing works with 2.x
```
pip install oftest conan==1.58.0
```

Now clone the repository:
```bash
git clone https://github.com/DLR-RY/FMU4FOAM.git
```

Before running `build-ECI4FOAM.sh`, remove `Allwmake` and `cd ..` from `build-ECI4FOAM.sh` at the end of the script. 
The reason is that the `Allwmake` will fail because of a couple of files that do not compile on new OF versions.

Now run: 
```bash
./build-ECI4FOAM.sh
```

Remove:
`$(inputBC)/coupledWallHeatFluxTemperature/coupledWallHeatFluxTemperatureFvPatchScalarField.C`
and
`$(output)/extForces/extForces.C`
from
`FMU4FOAM/ECI4FOAM/src/externalComm/Make/files`

Enter the `ECI4FOAM` and finally do the `Allwmake`

Get to the root folder of FMU4FOAM
```bash
cd ..
```

Run `./Allwmake`

Finally:
```bash
pip install fmu4foam
pip install OMSimulator
pip install pythonfmu
```

To test:
```bash
cd examples/heatedRoom
./Allrun
```




