# Testing for physics coupling

Tags: [![badge](https://img.shields.io/badge/ThermalHydraulics-onePhase-blue.svg)]() [![badge](https://img.shields.io/badge/Neutronics-diffusion-blue.svg)]() [![badge](https://img.shields.io/badge/ThermalMechanics-legacyThermoMechanics-blue.svg)]() [![badge](https://img.shields.io/badge/Multiphysics-looseCoupling-orange.svg)]()


## Description

This is a very simple 2-D case that you can use to test and explore physics coupling. It consists of a steady-state diffusion calculation and a transient point-kinetics calculation.

Each physics is simulated on a 2-D rectangular domain in the x-z plane:

- The N (neutronics) domain has one `cellZone` and spans from `x=0` to `x=0.75`.
- The TH (thermal-hydraulics) domain has one `cellZone` and spans from `x=0.25` to `x=0.5`.
- The TM (thermal-mechanics) domain has three `cellZones`:
  - one from `x=0.25` to `x=0.5`
  - one from `x=0.5` to `x=0.75`
  - one from `x=0.75` to `x=1`

```
|---|---|---|---|
| N   N   N |
    | TH|
    | TM  TM  TM|
```

This combination provides the TM domain with:

- a zone overlapped with both N and TH, where it can obtain the temperature from TH
- a region overlapped only with N, where it can obtain the `powerDensity` from N and compute its own temperature
- a region not overlapped with anything, where it computes `T` without sources

Similarly, the N domain will have:

- a region not overlapped with anything, where it uses default values of `XS`
- a region overlapped with both TM and TH, where it receives fields from both
- a region overlapped only with TM

In TH, coolant (properties ~ sodium) flows from the bottom of the domain to the top (along `z`).

In N, cross-sections have been arbitrarily chosen to obtain a critical system in two groups. They are used in an initial steady-state diffusion calculation to obtain the power shape used in the point-kinetics transient calculations.

The TM case features a fixed boundary on the right side.

If an initial diffusion calculation is not performed, the point-kinetics solver assumes a spatially uniform power.