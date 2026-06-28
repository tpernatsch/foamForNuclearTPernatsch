# 1D Boiling

Tags: [![badge](https://img.shields.io/badge/ThermalHydraulics-twoPhase-blue.svg)]()

## Description

This test case demonstrates the two-phase capabilities of GeN-Foam. It consists of a 1D channel with a pressure-driven flow of liquid sodium. A power source is turned on at time 0, which eventually leads to boiling. Because the flow is pressure driven, the system exhibits a flow excursion. After a certain time has elapsed (see `constant/phaseProperties.structureProperties.powerOffCriterionModel`), the power is turned off.

For further information, see the comments in [`constant/fluidRegion/phaseProperties`](constant/fluidRegion/phaseProperties), [`system/controlDict`](system/controlDict), and [`system/fluidRegion/fvSolution`](system/fluidRegion/fvSolution).