# Sub-solvers

Thermal-hydraulics calculations are performed by classes derived from *thermalHydraulicsModel.H* that contain specific sub-solvers:

- *onePhase* for single-phase calculations, using the formulation proposed in Refs. [@Radman2019ADesign] [@RADMAN2021111178] [@RADMAN2021111422]  (see [*onePhase.H*](../classes/thermalHydraulics/solvers/onePhase/onePhase.md))
- *twoPhase* for adjoint diffusion calculations, using the formulation proposed in Refs. [@Radman2019ADesign] [@RADMAN2021111178] [@RADMAN2021111422] (see [*twoPhase.H*](../classes/thermalHydraulics/solvers/twoPhase/twoPhase.md))

For the user, the derived classes translate into runtime selectable models. The specific sub-solver to be used in a simulation can be selected in the solvers dictionary like explained in the coupling section.


## OpenFOAM-based sub-solvers

The new GeN-Foam structure allows to include already developed OpenFOAM solver. This requires to transpose the solver application into a GeN-Foam solver class. *compressibleInterFoam* shows how to translate standard OpenFOAM solvers into this new format, taking into account all the dependencies.

We list below the imported OpenFOAM-based standard solvers available in GeN-Foam.

- *compressibleInterFoam* for two compressible, non-isothermal immiscible fluids using a VOF (volume of fluid) phase-fraction based interface capturing approach. This standard solver has been transposed from OpenFOAM ([link](https://www.openfoam.com/documentation/guides/latest/man/compressibleInterFoam.html)) into the new GeN-Foam solver structure (see [*src/classes/openFoamImportedSolvers/compressibleInterFoam*](../classes/openFoamImportedSolvers/compressibleInterFoam/compressibleInterFoam.md)).