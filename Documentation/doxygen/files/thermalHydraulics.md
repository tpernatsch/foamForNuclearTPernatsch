
# Thermal-hydraulics {#TH}

**Work in progress!!**

## Introduction

Both single- and two-phase simulations can be performed using GeN-Foam. All sub-solvers were developed for a coarse-mesh porous-medium treatment of complex structures such as core and heat exchanger, and for a standard RANS treatment of clear-fluid regions. The sub-solvers automatically switch from a porous-medium (coarse-mesh) treatment to a standard CFD (fine-mesh) treatment when the colume fraction of the sub-scale structures is set to zero. This allows for an implict coupling of porous-medium (sub-channel-like in 2D and 3D, or system-code-like) treatment of compelex structures (e.g., core and heat exchnagers) with a standard CFD treatment of clear-fluid regions (e.g., plena and pools).

A coarse-mesh porous-medium treatment of the core implies that the core is modeled without resolving the sub-scale structure (e.g., the fuel rods or the heat exchanger tubes). As a matter of fact, in principle and for consistency, the finest radial mesh chosen by a user should not finer than one cell per pin cell. A porous-medium formulation derives from a volume averaging of the Navier-Stokes equations. The volume averaging results in source terms that describe the interaction (drag and heat transfer) of the fluid with the sub-scale structure. In GeN-Foam, these source terms are modeled using user-selectable correlations for drag (e.g., correlations for the Darcy friction factor) and heat transfer (e.g., correlations for the Nusselt number). In this sense, a porous-medium model can be associated with a 3-D version of a system code.

With regards to the modelling of the sub-scale structures, GeN-Foam allows to model simultaneously in the same region both a "power model" and a "passive structure". Power models are used to model fo instance the nuclear fuel (based on a 1-D approximation), electrically heated rods, or a fixed temperature body (which can be used to approximate a heat exchanger). Passive structures are structures that passively heats up or cool down based on their own heat capacity, volumetric area, and heat tranfer with the coolant. This can be used to model structures like the assembly wrappers or the reflectors.

All thermal-hydraulics functionalities are handled by the class *thermalHydraulicsModel.H*, the derived classes for the various sub-solvers (see below), and a thermal-hydraulic library that can be found under */GeN-Foam/classes/thermalHydraulics/src*.

## Sub-solvers

Thermal-hydraulics calculations are performed by classes derived from *thermalHydraulicsModel.H* that contain specific sub-solvers:
* *onePhase* for single-phase calculations, using the formulation proposed in Refs. \cite Radman2019ADesign \cite RADMAN2021111178 \cite RADMAN2021111422  (see *onePhase.H*)
* *onePhaseLegacy* for single-phase calculations, using the formulation proposed in Ref. \cite FIORINA201524 (see *onePhaseLegacy.H*)
* *twoPhase* for adjoint diffusion calculations, using the formulation proposed in Refs. \cite Radman2019ADesign \cite RADMAN2021111178 \cite RADMAN2021111422 (see *twoPhase.H*)
For the user, the derived classes translate into runtime selectable models. The specific sub-solver  to be used in a simulation can be selected at runtime in the *phaseProperties* dictionary in *constant/fluidRegion/*, using the keyword *thermalHydraulicsType*. 


## Porous-medium properties

The various parameters to be used in a porous-medium simulation can be set using the *phaseProperties* dictionary. 
<div class="border-box" style='padding:0.1em; margin-left: 4em;  margin-right: 8em;  border: 1px solid gray; background-color:#f2f3fa; color:#05134a'>
<b>For the user: the *phaseProperties* dictionary</b>

The *phaseProperties* dictionary can be found in *constant/fluidRegion/*. It is a large dictionary that can be used to: chose the sub-solver to be used (one-phase, legacy one-phase or two-phase); set various properties of the phases (beside basic thermo-physical properties defined in the *thermophysicalProperties* dictionary); set the properties of the sub-scale structures (fuel pins, heat exchangers, etc) in the porous zones, including the possibility to assign a *powerModel* for power production (e.g., nuclear fuel, or constant power) and the *passiveProperties* of another sub-structure that interacts thermally with the fluid (for instance the wrappers in sodium fast reactors).  The name of the porous zones must coincide with that of the cellZones of the fluidRegion mesh. Anisotropic pressure drops can be set by using the keywords *transverseDragModel* (Blasius, GunterShaw, same) and *principalAxis*(localX, localY, localZ) in the sub-dictionary *dragModels.(nameOfPhase).structure.(nameOfCellZones)*. *principalAxis* sets the axis on which the nominal dragModel is used. *transverseDragModel* sets the model to be used on the two directions that are perpendicular to *principalAxis*. If *same* is chosen as *transverseDragModel*, the code will use the nominal model in all directions, but with the possibility of an anisotropic hydraulic diameter. The anisotropy of the hydraulic diameter can be set using the keyword *localDhAnisotrpy* and assign to it a vector of 3 scaling factors (one for each local directions). 
<br><br>One can find detailed, commented examples in the tutorials 
[3D_SmallESFR](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/3D_SmallESFR/rootCase/constant/fluidRegion/phaseProperties) (single phase) and
[1D_boiling](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/1D_boiling/constant/fluidRegion/phaseProperties) (two phases).
</p>
</div>
<br>

<div class="border-box" style='padding:0.1em; margin-left: 4em;  margin-right: 8em;  border: 1px solid gray; background-color:#f2f3fa; color:#05134a'>
<b>For the user: the *g* dictionary</b>


## Physical properties

The *g* dictionary can be found under *constant/fluidRegion/*. It is a standard OpenFOAM dictionary that allows specifying the gravitational acceleration.
</p>
</div>
<br>

<div class="border-box" style='padding:0.1em; margin-left: 4em;  margin-right: 8em;  border: 1px solid gray; background-color:#f2f3fa; color:#05134a'>
<b>For the user: the *thermophysicalProperties* dictionary</b>

The *thermophysicalProperties* dictionary can be found under *constant/fluidRegion/*. It is a standard OpenFOAM dictionary that allows defining the thermo-physical properties of the coolant. When performing two-phase flow analyses, two dictionaries must be employed named*thermophysicalProperties.(name of fluid)*. The name of the two fluids are defined in the *phaseProperties* dictionary.

One can find a detailed, commented example in hte tutorials 
[3D_SmallESFR](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/3D_SmallESFR/rootCase/constant/fluidRegion/thermophysicalProperties) (one-phase), 
[1D_boiling (liquid)](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/1D_boiling/constant/fluidRegion/thermophysicalProperties.liquid) (two-phase, liquid)
[1D_boiling (vapour)](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/1D_boiling/constant/fluidRegion/thermophysicalProperties.vapour) (two-phase, vapour)
</p>
</div>
<br>


## Turbulence properties

<div class="border-box" style='padding:0.1em; margin-left: 4em;  margin-right: 8em;  border: 1px solid gray; background-color:#f2f3fa; color:#05134a'>
<b>For the user: the *turbulenceProperties* dictionary</b>

The *turbulenceProperties* dictionary can be found under *constant/fluidRegion/*. It is a standard OpenFOAM dictionary that allows defining the turbulence model to be used.

One can find a detailed, commented example in the tutorial 
[3D_SmallESFR](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/3D_SmallESFR/rootCase/constant/fluidRegion/turbulenceProperties).
</p>
</div>
<br>

## Initial and boundary conditions

Initial and boundary conditions adopt the usual OpenFOAM logic for one- and two-phase solvers. A couple of things to be kept in mind:
<UL>
<LI> The pressure field we solve for is p_rgh (pressure minus the gravitational head)
<LI> When performing turbulent analyses, one need to add the fields *k*, *epsilon*, nut and *alphat*
</UL>

One thing that instead specific to GeN-Foam (except for the one-phase legacy sub-solver) and that one needs to keep in mind is that U (or u.(name of fluid)) are the real velocities, not the Darcy velocities. In a porous structure, they represent the actual velocity of the fluid, and not the velocity multiplied by the fluid fraction. For instance, U will increase when transiting from a high- to a low-porosity region.


OpenFOAM provides most of the boundary conditions one may need for thermal-hydraulics models. In addition, a few boundary conditions have been included in GeN-Foam and can be found in *GeN-Foam/classes/thermalHydraulics/src/boundaryConditions*. Information on the use of each boundary condition can be found in the header files (.H). 


## Discretization and solution

Details for discretization and solution of equations are handled in a standard OpenFOAM way, i.e., through the *fvSolution* and *fvSchemes* dictionaries in *constant/fluidRegion*. 

© All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 2021
