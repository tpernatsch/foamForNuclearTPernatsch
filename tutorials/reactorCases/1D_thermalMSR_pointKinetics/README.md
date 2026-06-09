# 1D thermal MSR Point-Kinetics

Tags: [![badge](https://img.shields.io/badge/ThermalHydraulics-onePhase-blue.svg)]() [![badge](https://img.shields.io/badge/Neutronics-pointKinetics-blue.svg)]() [![badge](https://img.shields.io/badge/Multiphysics-looseCoupling-orange.svg)]()


## Description

A 1D loop model of a thermal Molten Salt Reactor (MSR) with point-kinetics
neutronics. Unlike a fast MSR, the core contains a graphite moderator
structure, so heat is deposited in both the salt and the solid graphite.
The tutorial demonstrates:

- how to use the `lumpedParameterStructure` power model to distribute
  fission heat between the salt and the graphite
- how to set a power density in both fluid and structural regions
- liquid-fuel point-kinetics with delayed neutron precursor (DNP) transport:
  precursors are advected by the flowing salt around the full loop
- `initPrecursorsLiquidFuel`: automatic initialization of the precursor
  spatial distribution at the start of a transient
- pseudo-time-stepping precursor initialization (`pseudoTimeSteps`,
  `pseudoTimeStepSize`, `pseudoTimeRelaxation`): a backward-Euler pseudo-time
  loop that converges the precursor field to flow-equilibrium before physical
  time-stepping begins, giving a more accurate and robust starting point
- `fuelZones`: restriction of the DNP fission source to the fissile core
  zone only, preventing spurious precursor production in non-fissile loop
  sections (heat exchanger, pump, legs)

The reactor geometry consists of a 1D loop with five cell zones: `core`
(fissile), `HotLeg`, `ColdLeg`, `hx` (heat exchanger), and `pump`.


## Cases

The tutorial runs five sequential cases.

### steadyState

Eigenvalue calculation to establish the initial neutron flux and temperature
distribution. Runs from t = 0 to t = 100 s.

### transient

Point-kinetics transient with liquid-fuel delayed neutron precursor (DNP)
transport. Precursors are advected by the salt flow around the loop. At the
start of the transient (t = 100 s), `initPrecursorsLiquidFuel true` triggers
a steady-state solve to distribute the precursors consistently with the flow
field before time-stepping begins. Runs from t = 100 s to t = 400 s.

### transientEnd

One-second continuation of `transient` (t = 400 s to t = 401 s), used to
verify that the solution remains stable at the end of the transient.

### transient_substepping

Repeat of `transient` demonstrating the pseudo-time-stepping precursor
initialization feature. Two additions are made to `nuclearData`:

- `fuelZones ( core )` — restricts the DNP fission source to the fissile
  core zone only, preventing spurious precursor production in the loop
  sections.
- `pseudoTimeSteps 100; pseudoTimeStepSize 10.0; pseudoTimeRelaxation 1.0`
  — at the start of the transient, a backward-Euler pseudo-time loop runs
  100 steps of 10 s each (1000 s total pseudo-time >> 1/λ_min ≈ 80 s),
  driving the precursor spatial distribution to a well-converged
  flow-equilibrium state before the physical transient begins.

### transientEnd_substepping

One-second continuation of `transient_substepping` (t = 400 s to t = 401 s).
Cloned from `transient_substepping`, so `fuelZones` and pseudo-time-stepping
parameters are inherited. `initPrecursorsLiquidFuel` is set to `false` so
that the converged precursor state at t = 400 s is preserved rather than
re-initialized.


## Key nuclearData entries

| Entry | Default (`transient`) | Substepping cases |
|---|---|---|
| `initPrecursorsLiquidFuel` | `true` | `true` (`transient_substepping`), `false` (`transientEnd_substepping`) |
| `fuelZones` | absent (all cells) | `( core )` |
| `pseudoTimeSteps` | absent (steady-state fallback) | `100` |
| `pseudoTimeStepSize` | — | `10.0` |
| `pseudoTimeRelaxation` | — | `1.0` |
