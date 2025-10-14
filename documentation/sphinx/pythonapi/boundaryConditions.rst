.. _pythonapi_boundaryConditions:

.. module:: foamForNuclear.boundaryConditions

---------------------------------------------------------------
:mod:`foamForNuclear.boundaryConditions` -- Boundary Conditions
---------------------------------------------------------------

.. contents:: Table of Contents
    :local:


Basic Boundary Conditions
-------------------------

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template: myclassinherit.rst

    Patch
    Calculated
    FixedValue
    Mixed
    ZeroGradient


Constraint Boundary Conditions
------------------------------

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template: myclassinherit.rst

    Cyclic
    CyclicAMI
    Empty
    Wedge


Inlet Boundary Conditions
-------------------------

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template: myclassinherit.rst

    FixedFluxExtrapolatedPressure
    FixedFluxPressure
    FlowRateInletVelocity
    MappedFlowRate
    PressureInletOutletVelocity
    TotalPressure


Outlet Boundary Conditions
--------------------------

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template: myclassinherit.rst

    InletOutlet
    PressureInletOutletVelocity
    TotalPressure


Wall Boundary Conditions
------------------------

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template: myclassinherit.rst

    MappedField
    MappedFixedInternalValue
    MappedFixedValue
    MappedFlowRate
    MovingWallVelocity
    NoSlip
    Slip
    LumpedMassWallTemperature
    WallHeatTransfer


Coupled Boundary Conditions
---------------------------

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template: myclassinherit.rst

    FixedJump
    MappedField
    MappedFixedInternalValue
    MappedFixedValue
    MappedFlowRate
    MappedPatch
    MappedVelocityFluxFixedValue
    TemperatureCoupledBase


Generic Boundary Conditions
---------------------------

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template: myclassinherit.rst

    MappedField
    MappedFixedInternalValue
    MappedFixedValue
    PrghPressure
    PrghTotalHydrostaticPressure
    PrghTotalPressure
    Slip
    UniformFixedValue


GeN-Foam Boundary Conditions
----------------------------

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template: myclassinherit.rst

    AlbedoSP3
    EpsilonWallFunction
    FixedMassFlowRate
    KqRWallFunction
    NutkWallFunction


OFFBEAT Boundary Conditions
---------------------------

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template: myclassinherit.rst

    CoolantPressure
    FuelRodGap
    GapContact
    TractionDisplacement


FMU4FOAM Boundary Conditions
----------------------------

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template: myclassinherit.rst

    CoupledFlowRateInletVelocity
    CoupledFlowRateOutletVelocity
    CoupledUniformExternalValue


Custom Boundary Conditions
--------------------------

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template: myclassinherit.rst

    CustomPatch
