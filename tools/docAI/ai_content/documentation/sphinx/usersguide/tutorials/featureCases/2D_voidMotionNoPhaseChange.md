# 2D Void Motion No Phase Change

Tags: [![badge](https://img.shields.io/badge/ThermalHydraulics-twoPhase-blue.svg)]()

## Description

This is a very simple tutorial that demonstrates a two-phase case without mass transfer between phases. This is achieved by not inserting the `phaseChangeModel` subdictionary in the `phaseProperties` dictionary. For details on how to use the two-phase flow solver, please refer to the `1D_boiling` tutorial.

The main additional feature used in this tutorial, compared to `1D_boiling`, is the use of the `initialAlphas` subdictionary in the `vapourProperties` subdictionary of the `phaseProperties` dictionary. It is used to provide potentially different initial phase fractions for different `cellZones` (as an alternative to the use of `setFields`).