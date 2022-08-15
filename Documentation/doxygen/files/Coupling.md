
# Coupling and time stepping {#COUPLING}

**Work in progress!!**

## Introduction

GeN-Foam has been developed for the steady-state and transient analysis of reactors featuring pin-type, plate-type, or liquid fuel (viz., Molten Salt Reactors). It includes sub-solvers for neutronics, single- and two-phase thermal-hydraulics, and thermal-mechanics, with the choice of the physics to solve that can be made at runtime. Three different meshes are employed for neutronics, thermal-hydraulics, and thermal-mechanics. The selection of the physics to solve is done in *system/controlDict*


<div class="border-box" style='padding:0.1em; margin-left: 4em;  margin-right: 8em;  border: 1px solid gray; background-color:#f2f3fa; color:#05134a'>
<b>The *controlDict* dictionary</b>

The *controlDict* is an extended version of the one that is normally used in other OpenFOAM solvers. Compared to a standard OpenFOAM controlDict, it includes several keywords that allow to select:
<UL>
<LI> Which physics to solve via the keywords *solveFluidMechanics*, *solveEnergy*, *solveNeutronics*, solveThermalMechanics
<LI> The type of reactor via the keyword *liquidFuel*
<LI> The options for time stepping via the keywords *adjustTimeStep*, *maxDeltaT*, *maxCo* and *maxPowerVariation*, as well as *maxCoTwoPhase* and *marginToPhaseChange* for two-phase flow simulations.
</UL> An option for mesh manipulation called *removeBaffles*. This flag is not mandatory and allows to create a ghost thermal-hydraulics mesh without baffles.This ghost mesh allows for better mesh-to-mesh projections between different physics. WARNING: parallel execution not tested.

Fairly complete examples of *controlDict* for single-phase flow can be found in [2D_FFTF](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/blob/master/Tutorials/2D_FFTF/rootCase/system/controlDict) and [3D_SmallESFR](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/blob/master/Tutorials/3D_SmallESFR/rootCase/system/controlDict), while an explanation of the two-phase flow options can be found in [1D_boiling](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/blob/master/Tutorials/1D_boiling/system/controlDict).

</p>
</div>


The coupling between physics is achieved by projecting coupling variables from the mesh they are calculated, to the mesh they need to be used. The following figure shows the overall logic behind the coupling. 


\image html couplingFields.png width=300px


The thermal-hydraulics sub-solver operates on a certain computational domain \f$\Omega_{TH}\f$. Given an input volumetric power density \f$q\f$, the thermal-hydraulics sub-solver is tasked with predicting the resulting fluid temperature \f$T\f$, density \f$\rho_{fluid}\f$ and velocity \f$u\f$ fields, as well as relevant structure temperature fields \f$T_s\f$. For a two-phase treatment, the fields \f$\rho_{fluid}\f$, \f$T\f$, \f$u\f$ consist of mass-weighed mixture values. The velocity field \f$u\f$ is used for coupling only when simulating MSRs to advect the precursors. The field \f$q\f$ is the volumetric fuel power density and it can pertain either to a sub-scale structure (typically, the fuel rods) or the fluid itself (i.e.\ the  liquid fuel in MSRs), depending on the system under investigation. The symbol \f$T_s\f$ collectively denotes the temperature fields of the structures, which can range from the fuel and cladding of a nuclear fuel pin to control rod drivelines, wrappers, the diagrid, etc. This entirely depends on what the structure thermal models are supposed to represent in the cell zones where they have been defined.

The neutronics sub-solver operates on a computational domain \f$\Omega_N\f$ and is tasked with predicting the volumetric fuel power density \f$q\f$ for varying coupling fields. Not all of these fields are always used, depending on the selected type of neutronics treatment. In general terms, the diffusion, \f$S_N\f$, \f$SP_3\f$ treatments are capable of modeling reactivity feedbacks from: coolant temperature \f$T\f$ and density \f$\rho_{fluid}\f$, average fuel and cladding temperatures collectively denoted with \f$T_s\f$, fuel axial displacement and core radial displacement collectively denoted as \f$d\f$. As long as a parametrization of the macroscopic cross-sections against these quantities is provided, these feedbacks can be resolved. The feedback reactivities of the point-kinetics solver are described by standard feedback coefficients.

The thermal-mechanics sub-solver operates on a computational domain \f$\Omega_{TM}\f$ and is tasked with predicting an overall displacement field that can be used to deform all the computational domains \f$\Omega_{TH}\f$, \f$\Omega_{N}\f$, \f$\Omega_{TM}\f$.  The displacement field is decomposed into fuel axial displacement and core radial displacement fields collectively denoted as \f$d\f$, which are passed to the neutronics to model expansion-related feedbacks. 


<div class="border-box" style='padding:0.1em; margin-left: 4em;  margin-right: 8em;  border: 1px solid gray; background-color:#f2f3fa; color:#05134a'>
<b>Prioritization in coupling variables</b>

This is one of the most complex aspect of GeN-Foam.
...
</p>
</div>


The coupling between physics is obtained via fixed-point iterations and the time derivatives are solved based on a first-order implicit Euler scheme. The figure below reports  the overall coupling scheme. Within each time step, the outer iteration loop is used to couple the single-physics solution steps it encompasses. It can be performed a user-selected number of times or controlled by the convergence of the residuals of the slowest-converging physics.

\image html couplingAlgorithm.png width=500px


The parameters for the coupling are set in *system/fvSolution*.

<div class="border-box" style='padding:0.1em; margin-left: 4em;  margin-right: 8em;  border: 1px solid gray; background-color:#f2f3fa; color:#05134a'>
<b>The general *fvSolution* dictionary</b>

The general *fvSolution* dictionary is found under */system/*  and allows to specify parameters related to the coupling among physics, and in particular: the type of coupling (implicit or explicit); and the parameters that affect the tightness pof the implicit coupling.

A commented  *fvSolution* can be found in [3D_SmallESFR](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/blob/master/Tutorials/3D_SmallESFR/rootCase/system/fvSolution). 
</p>
</div>



© All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 2021
