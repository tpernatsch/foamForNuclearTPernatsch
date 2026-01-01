.. _userguide_thermalhydraulics_porousMedium:


the phaseProperties dictionary
------------------------------

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
   :maxdepth: 2

    FFdrag
