# 1D Heat Exchanger

Tags: [![badge](https://img.shields.io/badge/ThermalHydraulics-onePhase-blue.svg)]() [![badge](https://img.shields.io/badge/ThermalHydraulics-twoPhase-blue.svg)]()


## Description

These test cases demonstrate the use of the novel `heatExchanger` feature.
This feature thermally couples two mesh domains that exist within the same
mesh region, but whose cells are otherwise not connected to each other.

In particular, the test case consists of two parallel 1-D channels, each with
an inlet and an outlet. Axially, a heated `cellZone` exists only on one
channel. On each channel, two `cellZones` act as the heat exchanger primary
and secondary, respectively. This setup allows the heated fluid in the
heated channel to exchange heat with the cooler fluid flow in the unheated
channel.

The `onePhase` case demonstrates the feature for single-phase flow. The
`twoPhase` case demonstrates the feature for two-phase flow, in which the
vapour condenses in the `heatExchanger` `cellZone` and heats the fluid on
the other `heatExchanger` side.

`FunctionObjects` (defined at the end of the `controlDict` for both cases)
allow tracking quantities of interest at the end of each time step, such as
the liquid inlet and outlet temperature and the mass flow rate in both
channels. This enables verification of energy conservation.

The specific `heatExchanger` feature options related to its usage are
commented in the `constant/fluidRegion/phaseProperties` dictionary of each
case.