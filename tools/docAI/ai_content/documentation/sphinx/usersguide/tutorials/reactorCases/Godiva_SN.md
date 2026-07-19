# Godiva sphere using discrete ordinate SN

Tags: [![badge](https://img.shields.io/badge/Neutronics-SN-blue.svg)]()

## Description

This is a purely neutronics eigenvalue case that demonstrates how to use the discrete ordinate (SN) solver in GeN-Foam. It simulates the Godiva experiment, which consists of a small super-prompt-critical sphere of enriched uranium.

The file [`constant/neutroRegion/quadratureSet`](constant/neutroRegion/quadratureSet) contains a simple quadrature set with 4 directions per octant. A more complex (and more computationally demanding) quadrature set with 16 directions per octant can be found in [`constant/neutroRegion/extraQuadratureSets/quadratureSet16`](constant/neutroRegion/extraQuadratureSets/quadratureSet16). A simpler option, with 1 direction per octant, is available in [`constant/neutroRegion/extraQuadratureSets/quadratureSet1`](constant/neutroRegion/extraQuadratureSets/quadratureSet1).

The scattering anisotropy can be changed by modifying the `legendreMoments` flag in [`constant/neutroRegion/neutronicsProperties`](constant/neutroRegion/neutronicsProperties). In that case, you should ensure that the corresponding scattering matrices are provided in the same file. In the current tutorial, scattering matrices are provided up to the 5th moment. This makes the tutorial computationally intensive.

It is recommended to run the case using the [`Allrun_parallel`](./Allrun_parallel) bash script on a good computer. In principle, the SN solver could be used for time-dependent calculations. However, no acceleration techniques are currently implemented, so the solver is particularly slow. An [`Allclean`](./Allclean) script is provided to clean up the case.