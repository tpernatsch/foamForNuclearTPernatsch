# 2D-wedge ALFRED model with FMU balance of plant

Tags: [![badge](https://img.shields.io/badge/ThermalHydraulics-onePhase-blue.svg)]() [![badge](https://img.shields.io/badge/Neutronics-pointKinetics-blue.svg)]() [![badge](https://img.shields.io/badge/Multiphysics-looseCoupling-orange.svg)]() [![badge](https://img.shields.io/badge/Multicode-FMI-red.svg)]()

Author: Thomas Guilbaud, EPFL/Transmutex SA, 26/01/2023

---

The purpose of this case is to simulate a coupled GeN-Foam/FMU model of an entire power plant, from the core to the turbine. The reactor is the Lead Fast Reactor (LFR) ALFRED with a thermal power of 300 MWth.

![](./images/2D_LFRpowerPlant_LFRfull.png)

*Fig 1. Full plant modeling of the ALFRED Lead Fast Reactor using the FMI interface and Modelica.*


## GeN-Foam

Each case requires an FMU representing the balance of plant, modeled using Modelica. A first steady-state simulation is performed before each transient. The transient simulation relies on the point-kinetics sub-solver.

Before starting the simulation, make sure that your system meets the following requirements:
- OpenFOAM-v2412
- GeN-Foam/develop build with FMU4FOAM (`./Allwmake --fmi`)
- OpenModelica 1.21.0 (`omc --Version`)
- Python3.8 with PyFMI (2.10.0), Pandas (1.4.3), Numpy (1.23.1)


### FMU generation and coupling architecture

Recommendation: Save the provided FMUs if the `omc` command does not work.

Prerequisites:
- Install OpenModelica and the ThermoPower package

Several architectures have been tested with this model:
- GeN-Foam coupled with OpenModelica exported as an FMU
- GeN-Foam coupled with OpenModelica and Simulink exported as FMUs

Generate the FMU using the OpenModelica compiler:
```bash
omc FMUGenSecCir.mos # generate all FMUs
# or
omc FMUGenSecCirWithCtrl.mos # generate SecondaryCircuitWithControl.fmu from SecondaryCircuitWithControl.mo
# or
omc FMUGenSecCirWithoutCtrl.mos # generate SecondaryCircuitWithoutControl.fmu from SecondaryCircuitWithoutControl.mo
```

You can change the architecture in the [rootCase/system/controlDict](rootCase/system/controlDict) by modifying the variable `pyFileName`:
- `pyFileName secondaryCircuitModelica;` for GeN-Foam + OpenModelica
- `pyFileName secondaryCircuit2FMUs;` for GeN-Foam + OpenModelica + Simulink


### Load follow transient

The system relies on two PID controllers. The first controller acts on the grid frequency and controls the turbine valve admission to balance the mechanical and electrical power. The second PID controls the reactivity insertion or the beam modulation to improve the core power response.

![](images/2D_LFRpowerPlant_SecondaryCircuitWithControlLegendSource.png)

*Fig 2. Balance of plant model using Modelica used by GeN-Foam.*


### Run a simulation

This tutorial proposes two scenarios: one with load following using external reactivity, and one with an external source. For both simulations, a first GeN-Foam standalone simulation is performed to converge the neutronics and thermal-hydraulics meshes. This improves the stability of the subsequent coupled GeN-Foam/Modelica simulation (see [`./Allrun_GeN-FoamStandalone`](./Allrun_GeN-FoamStandalone)).

Then, each case performs a coupled GeN-Foam/Modelica steady-state simulation without feedback from the external reactivity or the external neutron source. Once both codes have converged in steady state, a final transient simulation is performed with a variation of electric power demand (load).

Run the load-follow transient:
```bash
./Allclean

# Optional, run in Allrun_loadFollowReactivity and Allrun_loadFollowSource
# ./Allrun_GeN-FoamStandalone

# In critical mode
./Allrun_loadFollowReactivity
# or in ADS mode
./Allrun_loadFollowSource
```

Analyse the data (example of use):
```bash
python3 Allplot_loadFollowReactivity.py transientLoadFollowReactivity/secondaryCircuit.csv
# or
python3 Allplot_loadFollowSource.py transientLoadFollowSource/secondaryCircuit.csv
# or
python3 Allplot_loadFollowReactivity.py transientLoadFollowSource/secondaryCircuit.csv transientLoadFollowReactivity/secondaryCircuit.csv
# or
python3 Allplot_FMU.py transientULOF/secondaryCircuit.csv
# or
python3 Allplot_FMUanimate.py transientULOF/secondaryCircuit.csv
# or
python3 Allplot_loadFollowReactivity.py transientBeamTrip/secondaryCircuit.csv
```

You are encouraged to explore the different architectures by changing the `pyFileName`.


### Results

![](./GeN-FoamAndModelica/results_allplot_loadFollowReactivity.png)

*Fig 3. Evolution of power plant parameters during a load-follow transient in critical configuration.*


## All FMUs !

Two special cases have been created in `GeN-FoamAs2FMU` and `GeN-FoamAs3FMU`, where GeN-Foam is embedded in an FMU. The architecture is more homogeneous.

Go to [GeN-FoamAs2FMU](GeN-FoamAs2FMU) or [GeN-FoamAs3FMU](GeN-FoamAs3FMU) and execute:
```bash
./Allrun
```

The results can be extracted as follows:
```bash
python3 Allplot.py
```

In the end, the results are the same as for GeN-Foam + OpenModelica (FMU) and GeN-Foam + OpenModelica (FMU) + Simulink (FMU).