# 2D Cavity Boussinesq

Tags: [![badge](https://img.shields.io/badge/ThermalHydraulics-onePhase-blue.svg)]()

## Description

This test case demonstrates the use of the Boussinesq feature. The Boussinesq
approximation is used for buoyancy-driven flows. It assumes that
temperature-induced density changes are relevant only in the calculation of
the buoyancy momentum source terms. In all other terms of the momentum
equation, the density is kept constant.

By selecting the `Boussinesq` `equationOfState` type in thermophysical
properties and providing the parameters `rho0`, `beta`, and `T0`, the code
employs the Boussinesq approximation as originally intended. Specifically,
`rho = rho0` everywhere except in the buoyancy term. The buoyancy term is
computed as a density gradient using
`rho = rho0*(1-beta*(T-T0))`. Although the Boussinesq `equationOfState` is a
standard OpenFOAM model, the actual implementation consistent with the
approximation (i.e., that uses `rho0` in place of `rho` where needed) is
specific to this solver.

In this case, the domain is a square with 10 cm sides. Two opposing walls are
maintained at 500 K and 1000 K. The fluid starts at rest at 500 K. The
cooled and heated wall patch normals are perpendicular to gravity, while
gravity acts in the positive X direction. The fluid has thermophysical
properties comparable to those of liquid sodium at 500 K.

A steady state is reached with the establishment of a circular liquid motion.
The fluid rises against gravity at the heated wall and descends at the
cooled wall.