# Compiling GeN-Foam {#COMPILE}

GeN-Foam is built as a standard OpenFOAM solver, based on the www.openfoam.com distribution. As such, one should first install the latest OpenFOAM version and then compile GeN-Foam by typing in a terminal in the GeN-Foam source code root folder:

```bash
# Clean GeN-Foam executable
./Allwclean
```

```bash
# Compile GeN-Foam
./Allwmake
# or
./Allwmake -j8
```

Please notice that a new version of OpenFOAM is released by ESI/OpenCFD twice a year. It may take a few weeks for the developers to update GeN-Foam to a new OpenFOAM release.


## Compiling GeN-Foam for FMU coupling

Compiling GeN-Foam for FMU coupling is explained on this [page](@ref FMU).
