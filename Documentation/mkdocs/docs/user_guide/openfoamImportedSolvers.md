# OpenFOAM-imported sub-solvers

The GeN-Foam structure allows to include already developed OpenFOAM solver. This requires to transpose the solver application into a GeN-Foam solver class. *compressibleInterFoam* shows how to translate standard OpenFOAM solvers into the solver format, taking into account all the dependencies.

We list below the imported OpenFOAM-based standard solvers available in GeN-Foam.

- *compressibleInterFoam* for two compressible, non-isothermal immiscible fluids using a VOF (volume of fluid) phase-fraction based interface capturing approach. This standard solver has been transposed from OpenFOAM ([link](https://www.openfoam.com/documentation/guides/latest/man/compressibleInterFoam.html)) into the GeN-Foam solver structure (see [*src/classes/openFoamImportedSolvers/compressibleInterFoam*](../classes/openFoamImportedSolvers/compressibleInterFoam/compressibleInterFoam.md)).