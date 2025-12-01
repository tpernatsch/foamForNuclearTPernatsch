# Godiva sphere using discrete ordinate SN

Tags: [![badge](https://img.shields.io/badge/Neutronics-SN-blue.svg)]()


## Description

This is a purely neutronics eigenvalue case displaying how to use the discrete ordinate (SN) solver of GeN-Foam. It simulates the Godiva experiment, constituted by a small super-prompt-critical sphere of enriched Uranium.

The SN solver is selected in the [`constant/neutroRegion/neutronicProperties`](constant/neutroRegion/neutronicProperties) dictionary. The file [`constant/neutroRegion/quadratureSet`](constant/neutroRegion/quadratureSet) contains a simple quadrature set with 4 directions per octant. A more complex (and more computationally requiring) quadrature set with 16 directions per octant can be found in [`constant/neutroRegion/quadratureSet16`](constant/neutroRegion/quadratureSet16). A simpler one, with 1 direction per octant, can be found in [`constant/neutroRegion/quadratureSet1`](constant/neutroRegion/quadratureSet1).

The scattering anisotropy can be changed by changing the `legendreMoments` flag in [`constant/neutroRegion/nuclearData`](constant/neutroRegion/nuclearData). Of course, one should make sure that the corresponding scattering matrices are provided in the same file. In the current tutorial, scattering matrices are provided till the 5th moment. This is a computationally intensive tutorial.

It is suggested to run it using the [`Allrun_parallel`](./Allrun_parallel) bash script on a good computer. In principle, the SN solver could be used for time-dependent calculations. However, no acceleration techniques are currently implemented, making the solver particularly slow. An [`Allclean`](./Allclean) script is provided to clean up the case.
