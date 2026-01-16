.. _userguide_thermalhydraulics_porousMedium:


The *phaseProperties* dictionary
--------------------------------

The *phaseProperties* dictionary is a specialized dictionary that must be used in solvers that feature one- or two-phase porous-medium capabilities: :ref:`onePhase <onePhase>`, :ref:`onePhaseLegacy <onePhaseLegacy>`, :ref:`twoPhase <twoPhase>`. As a reminder, porous-medium solvers can be seen as a 3-D generalization of system codes, where the fluid-structure interaction must be modeled through dedicated correlations for pressure drops and heat transfer, and where the behavior of the non-resolved structures require sub-scale models such as a 1-D model for nuclear fuel pins. In foamForNuclear, the basic assumption is that porous-medium correlations involve three "phases": two fluids (e.g., water and vapour) and one solid structure (e.g., the nuclear fuel). In one-phase simulations, only one fluid and one structure are modeled. When there is no structure, the model reverts back to standard free-flow CFD modeling. 

The *phaseProperties* dictionary is a complex dictionary that handles:
   - Properties of a homegenized structure in a porous-medium approximation (volume fraction, hydraulic diameter, thermal properties), including the possibility of power models, i.e., models that can be used to simulate the thermal behavior of power-producing structures such as nuclear fuel pins, electrically heated rods, pebbles, or more generic structtures that can be modeled through equivalent electric circuits.
   - Regime maps.
   - Physics models for drag (these can be traditional pressure drops correlations) and heat transfer (these can be traditional Nusselt number correlations). For two-phase flow simulations, drag and heat transfer models should also be specified for the different phase pairs (fluid 1 to to structure, fluid 2 to structure, fluid 1 to fluid 2). In addition, two-phase flow simulation require models for virtual mass coefficient, pair geometry models, phase change models
   - For two-phase flow simulations:
      - The name of the phases 
      - The propoerties of each phase (e.g., state of matter, dispersed diameter model, etc)
   - Some miscellaneous parameters such as ``residualKd``


Run-time selectable models
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. toctree::
   :maxdepth: 2

   powerModels
   FFdrag
   FSdrag
   multipliersDrag
   contactPartiotionModels
   dispersionModels
   fluidDimaterModels
   FFHeatTransferCoefficientModels
   interfacialAreaModels
   regimeMapModels
   virtualMassCoefficientModels
   phaseChangeModels


.. contactPartiotionModels
   dispersionModels
   fluidDimaterModels
   FFHeatTransferCoefficientModels
   FSHeatTransferCoefficientModels
   interfacialAreaModels
   phaseChangeModels ()
   regimeMapModels
   virtualMassCoefficientModels




Other models
~~~~~~~~~~~~~~~~~~~~~~~~~~

Pump
^^^^

Heat exchanger
^^^^^^^^^^^^^^

Power-off criterion
^^^^^^^^^^^^^^^^^^^^

.. The simplest structures in GeN-Foam are the passive structures. These passive
.. structures are simply modeled as a heat capacity and they can be used for
.. instance to model assembly wrappers or reflector structures. In essence, the
.. fluid will interact cell-by-cell with these passive structures: they will take
.. energy from the fluid if the temperature of the fluid is higher than that of the
.. surface of the structure, and vice versa. A more complex example of sub-scale
.. structure is given by the powerModels. An example of a power model is the
.. nuclearFuelPin, which can be used to model a standard pin-type fuel. This power
.. model is capable of getting the power density from neutronics, solving a 1-D
.. model for heat transfer in the fuel, and giving back to the fluid the
.. temperature at the surface of the cladding. The fluid will then be capable of
.. calculating the heat transfer with the fuel based on the cladding surface
.. temperature and the Nusselt number. Each cellZone can host one passive structure
.. and one powerModel.

.. Both single- and two-phase simulations can be performed using GeN-Foam. All
.. sub-solvers were developed for a coarse-mesh porous-medium treatment of complex
.. structures such as core and heat exchanger, and for a standard RANS treatment of
.. clear-fluid regions. The sub-solvers automatically switch from a porous-medium
.. (coarse-mesh) treatment to a standard CFD (fine-mesh) treatment when the volume
.. fraction of the sub-scale structures is set to zero. This allows for an implicit
.. coupling of porous-medium (sub-channel-like in 2D and 3D, or system-code-like)
.. treatment of complex structures (e.g., core and heat exchangers) with a standard
.. CFD treatment of clear-fluid regions (e.g., plena and pools).

.. A coarse-mesh porous-medium treatment of the core implies that the core is
.. modeled without resolving the sub-scale structure (e.g., the fuel rods or the
.. heat exchanger tubes). As a matter of fact, in principle and for consistency,
.. the finest radial mesh chosen by a user should not be finer than one cell per
.. pin cell. A porous-medium formulation derives from a volume averaging of the
.. Navier-Stokes equations. The volume averaging results in source terms that
.. describe the interaction (drag and heat transfer) of the fluid with the
.. sub-scale structure. In GeN-Foam, these source terms are modeled using
.. user-selectable correlations for drag (e.g., correlations for the Darcy friction
.. factor) and heat transfer (e.g., correlations for the Nusselt number). In this
.. sense, a porous-medium model can be associated with a 3-D version of a system
.. code.

.. With regards to the modeling of the sub-scale structures, GeN-Foam allows
.. modeling simultaneously in the same region both a "power model" and a "passive
.. structure". Power models are used to model for instance the nuclear fuel (based
.. on a 1-D approximation), electrically heated rods, or a fixed temperature body
.. (which can be used to approximate a heat exchanger). Passive structures are
.. structures that passively heat up or cool down based on their heat capacity,
.. volumetric area, and heat transfer with the coolant. This can be used to model
.. structures like assembly wrappers or reflectors.


