
# Coupling and time stepping {#COUPLING}

**Work in progress!!**

## Introduction

GeN-Foam has been developed for the steady-state and transient analysis of reactors featuring pin-type, plate-type, or liquid fuel (viz., Molten Salt Reactors). It includes sub-solvers for neutronics, single- and two-phase thermal-hydraulics, and thermal-mechanics, with the choice of the physics to solve that can be made at runtime. Three different meshes are employed for neutronics, thermal-hydraulics, and thermal-mechanics. The selection of the physics to solve is done in *system/controlDict*


<div class="border-box" style='padding:0.1em; margin-left: 4em;  margin-right: 8em;  border: 1px solid gray; background-color:#f2f3fa; color:#05134a'>
<b>The *controlDict* dictionary</b>

The *controlDict* is an extended version of the one that is normally used in other OpenFOAM solvers. Compared to a standard OpenFOAM controlDict, it includes several keywords that allow one to select:
<UL>
<LI> Which physics to solve via the keywords *solveFluidMechanics*, *solveEnergy*, *solveNeutronics*, solveThermalMechanics
<LI> The type of reactor via the keyword *liquidFuel*
<LI> The options for time stepping via the keywords *adjustTimeStep*, *maxDeltaT*, *maxCo* and *maxPowerVariation*, as well as *maxCoTwoPhase* and *marginToPhaseChange* for two-phase flow simulations.
<LI> An option for mesh manipulation called *removeBaffles*. This flag is not mandatory and allows to create a ghost thermal-hydraulics mesh without baffles.This ghost mesh allows for better mesh-to-mesh projections between different physics. WARNING: parallel execution not tested.
</UL>
Fairly complete examples of *controlDict* for single-phase flow can be found in [2D_FFTF](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/blob/master/Tutorials/2D_FFTF/rootCase/system/controlDict) and [3D_SmallESFR](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/blob/master/Tutorials/3D_SmallESFR/rootCase/system/controlDict), while an explanation of the two-phase flow options can be found in [1D_boiling](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/blob/master/Tutorials/1D_boiling/system/controlDict).
</div>

## Coupling logic

The coupling between physics is achieved by projecting coupling variables from the mesh they are calculated, to the mesh they need to be used. The following figure shows the overall logic behind the coupling. 

\image html CouplingGF.png width=500px

Given an input volumetric power density \f$Q\f$, the thermal-hydraulics sub-solver is tasked with predicting the resulting fluid temperature \f$T\f$, density \f$\rho\f$ and velocity \f$u\f$ fields, as well as relevant structure temperature fields \f$T_s\f$. For a two-phase treatment, the fields \f$\rho\f$, \f$T\f$, \f$u\f$ consist of mass-weighed mixture values. The velocity field \f$u\f$ is used for coupling only when simulating MSRs to advect the precursors. The field \f$Q\f$ is the volumetric fuel power density and it can pertain either to a sub-scale structure (typically, the fuel rods) or the fluid itself (i.e.\ the liquid fuel in MSRs), depending on the system under investigation. The symbol \f$T_s\f$ collectively denotes the temperature fields of the structures, which can range from the fuel and cladding of a nuclear fuel pin to control rod drivelines, wrappers, the diagrid, etc. This entirely depends on what the structure thermal models are supposed to represent in the cell zones where they have been defined.

The neutronics sub-solver is tasked with predicting the volumetric fuel power density \f$Q\f$ for varying coupling fields. Not all of these fields are always used, depending on the selected type of neutronics treatment. In general terms, the diffusion, \f$S_N\f$, \f$SP_3\f$ treatments are capable of modeling reactivity feedbacks from: coolant temperature \f$T\f$ and density \f$\rho\f$, average fuel and cladding temperatures collectively denoted with \f$T_s\f$, fuel axial displacement and core radial displacement collectively denoted as \f$d\f$, as well as the temperature predicted by the thermal-mechanics sub-solver (this is useful for heterogeneous or mixed treatments, see below). As long as a parametrization of the macroscopic cross-sections against these quantities is provided, these feedbacks can be resolved. The feedback reactivities of the point-kinetics solver are described by standard feedback coefficients.

The thermal-mechanics sub-solver is tasked with predicting temperatures in solid regions and an overall displacement field that can be used to deform the neutronics mesh.  The displacement field is decomposed into fuel axial displacement and core radial displacement fields collectively denoted as \f$d\f$, which are passed to the neutronics to model expansion-related feedbacks. With regards to temperature, it is forced to the temperature of the sub-scale structures of the thermal-hydraulics domain whenever there is an area of overlap between thermal-mechanics and thermal-hydraulics mesh. If there is no overlap, the temperature is calculated based on a simple heat diffusion equation with parameters specified in the *thermoMechanicalProperties* dictionary. If the area is overlapped with neutronics, it will take from there a volumetric power source.

Based on the above, one may guess that GeN-Foam can operate in two different modes, depending on how the temperature of structures is calculated and as shown in the figure below:
- Homogeneous / domain overlap: the thermal-hydraulics solver is responsible for calculating temperatures throughout the physical domain. This is the most typical case, where structures are assumed to be treated as sub-scale structures in a porous-medium treatment.
- Heterogeneous: the thermal-hydraulics solver and the thermal-mechanics solver are responsible for calculating temperatures in different parts of the domain. This could be used for instance when simulating a core with a porous-medium approach, and a large solid reflector using the thermo-mechanics solver. Or it could be used to simulate a fuel pin in a heterogeneous manner.

Of course, it is also possible to have hybrid approaches, where the temperature in one structure is calculated in part based on the temperatures predicted in the sub-scale structure by the thermal-hydraulics sub-solver, and partly by the thermo-mechanical solver itself (where there is not overlap with the thermo-hydraulics domain).

Tutorial *2D_fullCoupling* has been created to allow users to play around with the couplings and understand their logic. 

\image html HetHom.png width=500px

N.B.: There is no need for the thermal-hydraulics, thermal-mechanics, and neutronics domains to be the same. GeN-Foam will project fields in a clever way whenever there is no overlap between 2 domains.

N.B.2: It is possible not to solve for displacements in the thermo-mechanics sub-solver setting to *false* the keyword *solveDisplacement* in *thermoMechanicalProperties* .






## Prioritization in coupling variables

This is one of the most complex aspects of multi-physics solvers like GeN-Foam. First of all, one should understand that, in multi-physics simulations, the coupling variables \f$q\f$, \f$T\f$, \f$\rho\f$, \f$u\f$, \f$T_s\f$ and \f$d\f$ must exist both in their original form in the region (mesh) where they are created, and in their projected form in the region (mesh) where they are used. With that in mind, the question immediately arises about the origin of variables when multiple sub-solvers are activated.


In GeN-Foam, for nearly all the coupling fields, a sub-solver will look at the original form of the coupling variable. For instance, the neutronic sub-solver will read the temperature field of the fuel from the thermal-hydraulics sub-solver (i.e., from the *fluidRegion*). This is physically intuitive and allows  GeN-Foam to behave as expected in sequential runs. For instance, one will be able to run a thermal-hydraulic calculation first, and then start from there and automatically use in a neutronic calculation the resulting temperature and density fields in order to properly parametrize the cross-sections. 

However, a drawback exists for a user that wishes to execute a single-physics simulation, but to provide its own specific value to the coupling field. The standard behavior of OpenFOAM will force this user to provide this value in the region where the field is normally calculated. Following the example above, a user wishing to run a neutronic calculation with its own temperature and densities will have to generate a mesh for the *fluifRegion*, and provide the desired temperatures and density fields in the 0 (or other *startTime*) folder of that region.

While the mentioned drawback is expected to negatively impact a minority of GeN-Foam users, this is not the case for \f$q\f$ (i.e., the *powerDensity* field), for which it was decided to adopt a separate treatment. In this case, GeN-Foam will take the *powerDensity* field from the *fluiRegion* when neutronics is switched off, while the neutronic sub-solver will project its own *powerDensity* field to the *neutroRegion* when neutronics is switched on. This allows for stand-alone thermal-hydraulic simulations with a power density calculated e.g. using Serpent, without the need to create a mesh for the neutronics. At the same time, when neutronics is switched on, it will overwrite the *powerDensity* of the *fluidRegion* with its own calculated value. More information about how to set the initial powerDensity in GeN-Foam calculations is provided in the "Setting the *powerDensity*" box in [Thermal-hydraulics](@ref TH).

Of course, a consistent and unified approach could be obtained for all coupling variables by introducing additional flags to let the user decide where the fields must be taken from, and projected to, but for the moment it was decided that this would add unnecessary complexity to code and input deck.


## Coupling loop

The coupling between physics is obtained via fixed-point iterations and the time derivatives are solved based on a first-order implicit Euler scheme. The figure below reports the overall coupling scheme. Within each time step, the outer iteration loop is used to couple the single-physics solution steps it encompasses. It can be performed a user-selected number of times or controlled by the convergence of the residuals of the slowest-converging physics.

\image html couplingAlgorithm.png width=500px


The parameters for the coupling are set in *system/fvSolution*.

<div class="border-box" style='padding:0.1em; margin-left: 4em;  margin-right: 8em;  border: 1px solid gray; background-color:#f2f3fa; color:#05134a'>
<b>The general *fvSolution* dictionary</b>

The general *fvSolution* dictionary is found under */system/*  and allows to specify parameters related to the coupling among physics, and in particular: the type of coupling (implicit or explicit); and the parameters that affect the tightness of the implicit coupling.

A commented  *fvSolution* can be found in [3D_SmallESFR](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/blob/master/Tutorials/3D_SmallESFR/rootCase/system/fvSolution). 
</div>


