# Foam-For-Nuclear Python API

FoamForNuclear Python API to setup, run and post-process GeN-Foam/OFFBEAT simulation.
This API is tailored for nuclear application but could be used for other OpenFOAM-based solvers.


## How to install

```bash
pip3 install .
# or
python3 setup.py install --user
```

Execute the following command to check the installation.

```bash
pip3 show foamForNuclear
```


### Install with FMI

Need `FMU4FOAM>=1.1.0` and one of below list:
- `FMPy>=0.3.21`
- `OMSimulator>=2.1.1`
- `pyfmi>=2.11.0`


## How to uninstall

Execute the following command to uninstall the foamForNuclear package.

```bash
pip3 uninstall foamForNuclear
```


## To do (2025-05-08)

- [x] Add link to external polyMesh folders if no block mesh or polyMesh not present
- [x] Add link to external .unv files if no block mesh or polyMesh not present
- [p] Add link to external nuclearData files (see [`examples/reactorCases/3D_SmallESFR_inProgress/Allrun.py`](examples/reactorCases/3D_SmallESFR_inProgress/Allrun.py))
- [x] Add point-kinetics data in nuclearData
- [x] Add quadratureSet dict for neutronics
- [x] Add externalSource file
- [x] Generate nuclear data from Serpent
- [x] Generate nuclear data from OpenMC
- [ ] Finish phaseProperties little models
    - [x] Add pump
    - [x] Add heat exchanger
    - [x] Two phase models
    - [p] Add hydraulic parameters calculations based on pitch, element diameter, wire, ... just like Honeycomb
- [p] Finish offbeat solverDict (only the main structure)
- [p] Add functionObjects
- [ ] Add post-processing functions
    - [x] Plot over slice
    - [x] Plot mesh
    - [x] Extract over line
- [ ] Add FMI function object
    - [x] FMUSimulator
    - [ ] FMUController
    - [x] Add FMI coupling graph with Mermaid
    - [x] Add OpenFOAM export to FMU
    - [x] Add multi-FMU running with semi-implicit scheme
- [ ] Add field, file and dictionary readers
- [x] Test square lattice blockMesh
- [p] Add missing boundary conditions
- [x] Replace `setup.py` to `pyproject.toml`
- [ ] Create predefined Offbeat cases (1D, 2D, ...) with preset boundary conditions
- [ ] Create predefined mesh
    - [x] Pure wedge geometry
- [x] Add command to run the solver through the API (openmc.run()), able to run pre-processing commands (blockMesh, decomposePar, ...)
- [ ] Try to remove the point.id in BlockMesh.py
- [p] Handle multiple phase thermophysical properties
- [x] Finish PIMPLE control dict (see 1D_HX/onePhase tuto)
- [x] Add relaxationFactors dict in fvSolution
- [x] Add coupling and solving graph using Mermaid
    - [x] Add solving flowchart
- [x] Add compressibleInterFoam
- [ ] Add scalarTransportFoam
- [x] Shortcut to extract keff
- [x] Add a simple results animation in .gif format
- [x] Create lattice mesh
- [x] Add baffles to meshes
- [ ] Add 1D piping
    - [x] Create simple 1D pipe mesh and connections
    - [ ] Add porousBafflePressure cyclic BC to couple two pipes https://doc.openfoam.com/2306/tools/processing/boundary-conditions/rtm/derived/wall/porousBafflePressure/
    - [ ] Add a GUI for piping placement?
- [ ] Add default drag and heat transfer models per fluid
- [ ] Create a mesh exploded view
- [x] Create mesh
    - [x] Using BlockMesh
    - [p] Using a Voronoi graph to generate a full `polyMesh`
        - [ ] Simplify the creation of axial zones
    - [ ] Using the Gmsh Python API
- [ ]


## Logo

[link](https://patorjk.com/software/taag/#p=display&h=1&f=Slant&t=FFN%0A)


## Build Sphinx Documentation

```bash
# In the root folder of foam-for-nuclear
sphinx-build docs/source public
```
