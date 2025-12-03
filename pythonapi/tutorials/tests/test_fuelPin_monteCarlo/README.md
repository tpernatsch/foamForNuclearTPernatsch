# Fuel pin and MGXS flow

Tags: [![badge](https://img.shields.io/badge/Neutronics-diffusion-blue.svg)]() [![badge](https://img.shields.io/badge/Neutronics-SP3-blue.svg)]()


## Description

This is an example of usage of MGXS generation and use in foamForNuclear using
the Python API.

The folder structure is as follow:

- [Serpent](./Serpent/) folder containing the Serpent model,
- [OpenMC](./OpenMC/) folder containing the OpenMC model,
- [GeN-Foam](./GeN-Foam/) folder containing the GeN-Foam model.


## GeN-Foam

The `GeN-Foam` folder contains 2 scripts:

- [`Allrun.py`](./GeN-Foam/Allrun.py) to run a neutronics calculation with one of the MGXS data set (default from Serpent).
- [`Allrun_datasets.py`](./GeN-Foam/Allrun_datasets.py) to run with multiple solvers and neutronics nuclearData sets.


### How to run

```bash
python3 Allclean.py

python3 Allrun.py
# or
python3 Allrun_datasets.py
```
