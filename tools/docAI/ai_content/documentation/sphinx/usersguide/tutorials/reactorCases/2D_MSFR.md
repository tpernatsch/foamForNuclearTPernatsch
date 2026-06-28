# 2D Molten Salt Fast Reactor

Tags: [![badge](https://img.shields.io/badge/ThermalHydraulics-onePhase-blue.svg)]() [![badge](https://img.shields.io/badge/Neutronics-diffusion-blue.svg)]() [![badge](https://img.shields.io/badge/Multiphysics-picardLoop-orange.svg)]()


## Description

`2D_MSFR` is a 2-D r-z model of a Molten Salt Fast Reactor. It solves for neutronics and thermal-hydraulics.

The `Allrun` bash script can be used to run the tutorial. The `Allclean` bash script can be used to clean it up. The script first runs a steady-state case with fluid dynamics only. Starting from the results of this simulation, a second steady-state case is launched to solve the neutronics and energy equations. Finally, a simple transient calculation is run. No reactivity is inserted in the transient, so the power remains constant for 10 seconds.

It is also possible to run one or more cases explicitly by passing their names as arguments to the script. This can be useful for hands-on sessions and learning. For example, the script can be used to run only the first case, and students can then prepare and launch the coupled and transient simulations themselves. Examples:

```bash
./Allrun steadyStateTH
./Allrun steadyStateEN transient
```

Valid case names are: `steadyStateTH`, `steadyStateEN`, `transient`. If no arguments are given, all three cases are executed in sequence. Invalid arguments will result in an error message.

Any modification to the initial conditions of the transient case will trigger an actual transient. For example, modifying `keff` in the [`reactorState`](./rootCase/0/uniform/reactorState) dictionary will trigger a reactivity-initiated transient. A more realistic transient can be initiated by modifying the heat transfer in the heat exchanger in the `phaseProperties` dictionary. The case is similar to the one presented in Ref. [1]. Please note that, to reduce computing time, the fluid-dynamics equations are not solved in the second steady-state and in the transient simulation. Also note that an upwind scheme is employed for the divergence term in the diffusion equations (in [`system/neutroRegion/fvSchemes`](./rootCase/system/neutroRegion/fvSchemes)), which is necessary to achieve convergence.

Lastly, a considerably finer mesh is provided under `constant/*/polyMeshFiner` directories.


## Run

```bash
./Allclean

./Allrun
```


## Gallery

<img src="images/2D_MSFR_Q.png" alt="Power density" width="450"/>

*Fig 1. Power density distribution.*


<img src="images/2D_MSFR_T.png" alt="Fluid temperature distribution" width="450"/>

*Fig 2. Fluid temperature distribution.*


<img src="images/2D_MSFR_Prec6.png" alt="Power density" width="450"/>

*Fig 3. 6th precursor group distribution.*