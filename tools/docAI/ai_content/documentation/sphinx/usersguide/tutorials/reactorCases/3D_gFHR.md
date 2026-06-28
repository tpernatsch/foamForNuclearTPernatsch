# 3D gFHR Pebble Bed Reactor

Tags: [![badge](https://img.shields.io/badge/ThermalHydraulics-onePhase-blue.svg)]()

This tutorial is mainly prepared based on the work of Yves Robert and Ludovic Jantzen at UC Berkeley. For further details, see: Robert, Y., et al., 2023. *"IMPACT OF THERMAL COUPLING ON A PEBBLE BED REACTOR EQUILIBRIUM FROM HYPER-FIDELITY DEPLETION"*, Proceedings of the 2023 30th International Conference on Nuclear Engineering, ICONE30, May 21–26, 2023, Kyoto, Japan.

## Description

This is a simple 3D model of the Kairos gFHR. It shows how to use the *nuclearSteadyStatePebble* and *lumpedNuclearStructure* power models to describe pebbles and TRISO.

The *nuclearSteadyStatePebble* model is a steady-state model purpose-built for pebble bed reactors. It allows you to include several details of the pebble/TRISO geometry and composition in an intuitive way.

The *lumpedNuclearStructure* model is a general model that lets you set up an electrical equivalent made of resistances and capacitances arranged in series. It can be used for transient calculations, but it requires you to translate your problem into an electrical equivalent. For this tutorial, such a derivation is provided in the *derivations* document, with the calculations performed in the *lumped_structure* Python script.

The `Allrun` file simply runs two different steady-state cases using the two models. You can check that the two models are equivalent, for instance, by plotting the maximum TRISO temperatures in ParaView.

## How to run

To remove previous results, if any:
```bash
./Allclean
```

To automatically run the two steady-state cases:
```bash
./Allrun
```