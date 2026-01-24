.. _userguide_thermalhydraulics_porousMedium:


The *phaseProperties* dictionary
--------------------------------

The *phaseProperties* dictionary is a specialized dictionary that must be used in solvers that feature one- or two-phase porous-medium capabilities: :ref:`onePhase <onePhase>`, :ref:`onePhaseLegacy <onePhaseLegacy>`, :ref:`twoPhase <twoPhase>`. A porous-medium formulation derives from a volume averaging of the Navier-Stokes equations. The volume averaging results in source terms that
describe the interaction (drag and heat transfer) of the fluid with the
sub-scale structure. In foamForNuclear, these source terms are modeled using
user-selectable correlations for drag (e.g., correlations for the Darcy friction
factor) and heat transfer (e.g., correlations for the Nusselt number). In this
sense, a porous-medium model can be associated with a 3-D version of a system
code. The non-resolved structures require sub-scale models such as a 1-D model for nuclear fuel pins. In foamForNuclear, the  assumption is that porous-medium simulations involve up to three "phases": two fluids (e.g., water and vapour) and one solid structure (e.g., the nuclear fuel). In one-phase simulations, only one fluid and one structure are modeled. In cell zones where there is no structure, the model reverts back to standard OpenFOAM free-flow CFD modeling, neabling hybrid simulations where certain components (e.g., a core) are modeled based on a porous-medium approximation, and others (e.g., the plena) are modeled using standard CFD.

The *phaseProperties* dictionary is a complex dictionary that handles:
   - Properties of a homegenized structure in a porous-medium approximation (volume fraction, hydraulic diameter, thermal properties), including the possibility of power models, i.e., models that can be used to simulate the thermal behavior of power-producing structures such as nuclear fuel pins, electrically heated rods, pebbles, or more generic structtures that can be modeled through equivalent electric circuits.
   - Regime maps.
   - Physics models for drag (these can be traditional pressure drops correlations) and heat transfer (these can be traditional Nusselt number correlations). For two-phase flow simulations, drag and heat transfer models should also be specified for the different phase pairs (fluid 1 to to structure, fluid 2 to structure, fluid 1 to fluid 2). In addition, two-phase flow simulation require models for virtual mass coefficient, pair geometry models, phase change models
   - For two-phase flow simulations:
      - The name of the phases
      - The properties of each phase (e.g., state of matter, dispersed diameter model, etc)

Detailed, commented examples are provided in the tutorials
`3D_SmallESFR <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/3D_SmallESFR/extendedThermoMechanics/constant/fluidRegion/phaseProperties>`_ (single phase) and
`1D_boiling <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/featureCases/1D_boiling/constant/fluidRegion/phaseProperties>`_
(two phases). In addition, an example of how to use a two-dimensional
flow-regime map can be found in
`1D_PSBT_SC <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/featureCases/1D_PSBT_SC/Phase_Ex1_12223/constant/fluidRegion/phaseProperties>`_.


Sub-dictionaries
~~~~~~~~~~~~~~~~

.. toctree::
   :maxdepth: 1

   structureProperties/index
   regimeMapModels/index
   physicsModels/index
   physicsModels/indexTwoPhase
   twoPhaseSpecific/index
   ../../../../cppapi/generated/porousMediaModels/phaseModels/structureModels/pump/pump.rst
   ../../../../cppapi/generated/porousMediaModels/phaseModels/structureModels/heatExchanger/heatExchanger.rst










.. Both single- and two-phase simulations can be performed using GeN-Foam. All
.. sub-solvers were developed for a coarse-mesh porous-medium treatment of complex
.. structures such as core and heat exchanger, and for a standard RANS treatment of
.. clear-fluid regions. The sub-solvers automatically switch from a porous-medium
.. (coarse-mesh) treatment to a standard CFD (fine-mesh) treatment when the volume
.. fraction of the sub-scale structures is set to zero. This allows for an implicit
.. coupling of porous-medium (sub-channel-like in 2D and 3D, or system-code-like)
.. treatment of complex structures (e.g., core and heat exchangers) with a standard
.. CFD treatment of clear-fluid regions (e.g., plena and pools).
