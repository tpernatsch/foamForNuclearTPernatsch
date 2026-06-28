# 1D MSR Point-Kinetics

Tags: [![badge](https://img.shields.io/badge/ThermalHydraulics-onePhase-blue.svg)]() [![badge](https://img.shields.io/badge/Neutronics-pointKinetics-blue.svg)]() [![badge](https://img.shields.io/badge/Multiphysics-looseCoupling-orange.svg)]()


## Description

This tutorial shows how to use the point kinetics module of GeN-Foam for
MSRs. It is a simple 1-D case consisting of a core, hot leg, pump, heat
exchanger, and cold leg. The geometry is one-dimensional. Salt
recirculation is simulated using a cyclic boundary condition between the
top and bottom boundaries.

Three simulations are performed:

- First, energy and fluid dynamics are coupled to obtain a steady state at
  the desired power level.
- Second, energy, fluid dynamics, and point kinetics are coupled to
  simulate a loss-of-flow.
- Third, the simulation is run to allow GeN-Foam to recalculate the
  reactivity loss caused by the recirculation of delayed neutron
  precursors. This can be used to verify the results against analytical
  solutions.

Please note that the power for the steady state is imposed via the
`powerDensity` field specified in the `fluidProperties` sub-dictionary in
the `constant/fluidRegion/phaseProperties` dictionary. It models a
constant power density in the core and zero power density in the other
regions. This power is then automatically updated by the point kinetics
solver during the transient.

Please also note that a correct evaluation of the reactivity worth of
delayed neutron precursors in MSRs would require knowledge of the adjoint
flux. In GeN-Foam, the adjoint flux is approximated by the `oneGroupFlux`.
When fluxes are not calculated via a diffusion calculation, you must
manually provide the `oneGroupFlux` in `0/neutroRegion`. In this tutorial,
the `oneGroupFlux` has been set to 1 in the core and zero elsewhere. This
is done via the `initialOneGroupFluxByZone` keyword in `nuclearData`.

A few Python files are provided to plot essential results (for instance:
`python plotPKPower_Temp.py ./transient/log.GeN-Foam`). In addition, a
`.m` file (that can be run using Octave) is provided to calculate expected
results at the end of the transient.

The tutorial was prepared by Arnaldo Mattioli (Politecnico di Milano) and
revised by Carlo Fiorina (EPFL) and Stefan Radman (EPFL).