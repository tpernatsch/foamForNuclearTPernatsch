
# Thermal-mechanics {#TM}

**Work in progress!!**

## Introduction

The thermal-mechanics solver of GeN-Foam is a simple linear elasticity solver that can be used to evaluate thermal deformations in a core. Temperatures are projected from the thermal-hydraulic solver (temperatures of fuels and structures) and the deformation field is used to deform the mesh for neutronics and thermal-hydraulics. In particular, the radial deformation of structures and the axial deformation of fuel are employed. to deform the neutronics mesh. 


## Models



## Various properties


## Initial and boundary conditions

Besides the standard ones available in OpenFOAM, GeN-Foam includes a *tractionDisplacement* boundary condition that allows to set a pressure or a traction on a boundary. Work is ongoing to adapt to GeN-Foam the contact boundary condition of the OFFBEAT fuel behavior solver \cite SCOLARO2020110416. 


## Discretization and solution

Details for discretization and solution of equations are handled in a standard OpenFOAM way, i.e., through the *fvSolution* and *fvSchemes* dictionaries in *constant/thermoMechanicalRegion*. 

© All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 2021
