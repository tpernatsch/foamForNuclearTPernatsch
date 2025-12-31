.. _userguide_thermalhydraulics:


Thermal-hydraulics
==================

Introduction
------------

Both single- and two-phase simulations can be performed using GeN-Foam. All
sub-solvers were developed for a coarse-mesh porous-medium treatment of complex
structures such as core and heat exchanger, and for a standard RANS treatment of
clear-fluid regions. The sub-solvers automatically switch from a porous-medium
(coarse-mesh) treatment to a standard CFD (fine-mesh) treatment when the volume
fraction of the sub-scale structures is set to zero. This allows for an implicit
coupling of porous-medium (sub-channel-like in 2D and 3D, or system-code-like)
treatment of complex structures (e.g., core and heat exchangers) with a standard
CFD treatment of clear-fluid regions (e.g., plena and pools).

A coarse-mesh porous-medium treatment of the core implies that the core is
modeled without resolving the sub-scale structure (e.g., the fuel rods or the
heat exchanger tubes). As a matter of fact, in principle and for consistency,
the finest radial mesh chosen by a user should not be finer than one cell per
pin cell. A porous-medium formulation derives from a volume averaging of the
Navier-Stokes equations. The volume averaging results in source terms that
describe the interaction (drag and heat transfer) of the fluid with the
sub-scale structure. In GeN-Foam, these source terms are modeled using
user-selectable correlations for drag (e.g., correlations for the Darcy friction
factor) and heat transfer (e.g., correlations for the Nusselt number). In this
sense, a porous-medium model can be associated with a 3-D version of a system
code.

With regards to the modeling of the sub-scale structures, GeN-Foam allows
modeling simultaneously in the same region both a "power model" and a "passive
structure". Power models are used to model for instance the nuclear fuel (based
on a 1-D approximation), electrically heated rods, or a fixed temperature body
(which can be used to approximate a heat exchanger). Passive structures are
structures that passively heat up or cool down based on their heat capacity,
volumetric area, and heat transfer with the coolant. This can be used to model
structures like assembly wrappers or reflectors.

All thermal-hydraulics functionalities are handled by the class
:ref:`thermalHydraulicsModel.H <thermalHydraulicsModel>`, the derived classes for the various sub-solvers
(see below), and a thermal-hydraulic library that can be found under
`foamForNuclear/src/porousMediaModels <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/src/porousMediaModels>`_.


The porous-medium approach in GeN-Foam
--------------------------------------

GeN-Foam was born for safety analyses and, to reduce computational footprint,
its base approach is to model for instance the core as a porous medium. In a
porous-medium approach, the fuel and other structures (e.g., the assembly
wrappers) are modeled using sub-scale models. This means that, in each cell, we
have the fluid and one or two lumped models for the sub-scale structures. The
simplest structures in GeN-Foam are the passive structures. These passive
structures are simply modeled as a heat capacity and they can be used for
instance to model assembly wrappers or reflector structures. In essence, the
fluid will interact cell-by-cell with these passive structures: they will take
energy from the fluid if the temperature of the fluid is higher than that of the
surface of the structure, and vice versa. A more complex example of sub-scale
structure is given by the powerModels. An example of a power model is the
nuclearFuelPin, which can be used to model a standard pin-type fuel. This power
model is capable of getting the power density from neutronics, solving a 1-D
model for heat transfer in the fuel, and giving back to the fluid the
temperature at the surface of the cladding. The fluid will then be capable of
calculating the heat transfer with the fuel based on the cladding surface
temperature and the Nusselt number. Each cellZone can host one passive structure
and one powerModel.

.. toctree::
   :numbered:
   :maxdepth: 3

    porousMedium/index

