# FoamForNuclear Python API

The foamForNuclear Python API is an interface designed to simplify the use of the OpenFOAM-based solvers GeN-Foam and OFFBEAT for multi-physics reactor and fuel performance simulations. More information in the [online documentation](https://foamfornuclear.gitlab.io/foamForNuclear/).


## Getting started

```bash
# Clone the repo
git clone --recursive https://gitlab.com/foamForNuclear/foamForNuclear.git

# Install the Python API
cd foamForNuclear/pythonapi
pip3 install -e .

# Or
cd foamForNuclear/
./Allwmake -j<N> --api
```


## Tutorials

This API comes with several [tutorials](./tutorials/) that reproduces the [standard tutorial](./../tutorials/) cases of the foamForNuclear project.
