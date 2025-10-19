# foamForNuclear (FFN)

**foamForNuclear (FFN)** is a **general-purpose, OpenFOAM®-based multi-physics platform** for the analysis and design of nuclear systems. It originates from the **merging and extension** of two OpenFOAM-based nuclear projects: **GeN-Foam**, dedicated to multi-physics reactor simulation, and **OFFBEAT**, focused on advanced nuclear fuel performance modeling.

FFN provides a **modular and extensible framework** for simulating a wide range of physics, from core neutronics and thermal-hydraulics to advanced thermomechanics and fuel behavior. Each physics capability is implemented as an independent **module**, solving a specific set of equations that can run standalone or be coupled with others in a multi-physics environment.


## Physics Modules

FFN integrates several physics modules developed through extensive research in nuclear system modeling:

- **Neutronics:** point kinetics, diffusion, adjoint diffusion, SP3, and discrete ordinates (SN) (steady or transient);
- **Thermal-hydraulics:** one-phase RANS-CFD and porous-medium models, and a two-phase porous-medium Euler–Euler model for sodium and water;
- **Solid temperature models:** for sub-scale structures, including 1-D fuel, fixed temperature/power, heated rods, fuel pebbles, and lumped-parameter models;
- **Thermomechanics:** linear/nonlinear elasticity, plasticity, creep, and temperature-dependent material properties;
- **Fuel behavior:** densification, swelling, fission gas release, creep, irradiation growth, non-conformal gap heat transfer, and burnup-dependent material properties.

The **modular structure** of FFN allows users familiar with the OpenFOAM® API to easily add new solvers to the list of available physics. As an example, a **Volume of Fluid** solver from OpenFOAM® is already integrated into FFN, enabling free-surface flow simulations.


## Applications within FFN

The FFN platform currently features two primary solver applications that leverage its physics modules: GeN-Foam and OFFBEAT. GeN-Foam offers flexible integration, allowing users to combine any set of physics modules—including OFFBEAT itself—for comprehensive simulations. For convenience, a standalone version of OFFBEAT is also available. Additionally, users familiar with OpenFOAM can easily develop custom solvers tailored to specific needs. This is particularly useful for specialized applications that require only one or two physics modules, where the full capabilities of GeN-Foam may be unnecessary.


### GeN-Foam

GeN-Foam is a **multi-physics application** enabling the coupling of multiple modules within a single simulation. GeN-Foam provides the infrastructure for:

- Defining **independent physics regions** (e.g., neutronics, thermal-hydraulics, thermomechanics, fuel behavior);
- Choosing **arbitrary coupling types** (volume, surface or interpolation-based for any field);
- Customizing **time-loops**, allowing both loosely and tightly coupling.

This architecture supports multi-scale analyses combining coarse and detailed models in a unified simulation.


### OFFBEAT

While the fuel behavior module can be integrated into GeN-Foam, **OFFBEAT** remains a **standalone application** within FFN, dedicated to **single-mesh fuel performance** simulations. It provides high-fidelity thermo-mechanical and material evolution modeling for individual fuel pins, rods, or pebbles.


## Documentation

Resources for users and developers include:

- **User Guide and Theory Manual:** [![foamForNuclear User's Guide](https://img.shields.io/badge/foamForNuclear-User_Guide-blue?logo=sphinx)](https://foamfornuclear.gitlab.io/foamForNuclear/index.html)
- **Online Doxygen API:** [![foamForNuclear Doxygen](https://img.shields.io/badge/foamForNuclear-Doxygen-blue?logo=doxygen)](https://foamfornuclear.gitlab.io/foamForNuclear/doxygen/index.html)
- **Python API Documentation:** [![foamForNuclear Doxygen](https://img.shields.io/badge/foamForNuclear-Python_API-blue?logo=python)](https://foamfornuclear.gitlab.io/foamForNuclear/pythonapi/index.html)
- **Introductory Lectures** ([`documentation/usefulDocumentsAndPresentations/`](./documentation/usefulDocumentsAndPresentations/))
- **Tutorial Cases** for each physics module and coupling type ([`tutorials`](./tutorials/))

Users are also encouraged to make use of the typical OpenFOAM learning strategies:

- the high-level C++-based object-oriented language of OpenFOAM, which normally allows understanding the logic of a solver easily;
- the comments that are typically available in the source code and, in particular, in the header files of each class;
- the support of the community.


## OpenFOAM Version

FFN is based on the **OpenFOAM® (ESI/OpenCFD)** distribution, currently **v2506**, available at [www.openfoam.com](https://www.openfoam.com). The platform is regularly updated to maintain compatibility with new releases.

## Getting started

First, [install OpenFOAM](https://www.openfoam.com/). If installing from source, [follow the instructions here](https://develop.openfoam.com/Development/openfoam/-/blob/master/doc/Build.md). Check for the appropriate version specified above!

Make sure OpenFOAM has been installed correctly and that the environment is correctly sourced. For instance, when installing from source, and if installing under the ~/openfoam folder:
```bash
source ~/openfoam/OpenFOAM-v2506/etc/bashrc
```

Once you are sure the appropriate OpenFOAM version is correctly installed, and the OpenFOAM environment is sourced:

```bash

# Clone the repo
git clone --recursive https://gitlab.com/foamForNuclear/foamForNuclear.git

# Compile the foamForNuclear project and build the foamForNuclear Python API
cd foamForNuclear
./Allwmake -j4 --api
```
j4 is telling your compiler to use 4 cores. You can use as many as your system allows.


## Copyright

© Contributions are individually acknowledged in the header files.


## Gallery

*Modeling of the European Sodium Fast Reactor: Boiling in a windowed assembly and core flowering*

<div style="text-align:center;">
  <img src="./documentation/media/assmblyWindows.png" height="200">
  <img src="./documentation/media/coreFlowering.png" height="200">
</div>


*Full-plant modeling of the ALFRED Lead Fast Reactor using the FMI interface and Modelica*

<div style="text-align:center;">
  <img src="./documentation/media/LFRfull.png" width="1000">
</div>


*Modeling of Molten Salt Reactors: MSRE and MSFR*

<div style="text-align:center;">
  <img src="./documentation/media/MSRE.png" height="300">
  <img src="./documentation/media/precTot3D.png" height="300">
  <img src="./documentation/media/precTot3Dside.png" height="300">
</div>


*Modeling of FFTF: 2-D primary circuit thermal-hydraulics and core fluxes*

<div style="text-align:center;">
  <img src="./documentation/media/FFTF.png" width="800">
</div>
