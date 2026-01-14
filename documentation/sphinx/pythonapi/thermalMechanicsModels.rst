.. _pythonapi_thermalMechanicsModels:

.. module:: foamForNuclear.thermomechanicalMaterial

-----------------------------------------------------------------------------
:mod:`foamForNuclear.thermomechanicalMaterial` -- Thermo-mechanical Materials
-----------------------------------------------------------------------------

.. contents:: Table of Contents
    :local:


Base and abstract thermo-mechanical class
=========================================

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template: myclassinherit.rst

    BaseThermomechanicalMaterial
    ThermomechanicalPropertyModel
    ConductivityModel
    DensityModel
    EmissivityModel
    HeatCapacityModel
    PoissonRatioModel
    ThermalExpansionModel
    YoungModulusModel
    DensificationModel
    SwellingModel
    PhaseTransitionModel
    RelocationModel
    FailureModel
    PoreVelocityModel
    YieldStressModel
    HardningYieldStressModel
    CreepModel
    LimbackCreepModel
    RheologyConstitutiveLaw
    ElasticityRheologyModel
    MisesPlasticityRheologyModel
    MisesPlasticCreepRheologyModel
    MaterialModel
    FuelMaterialModel
    ConstantMechanicalLaw


Constant Material
=================

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template: myclassinherit.rst

    ConstantMaterial


UO2
===

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template: myclassinherit.rst

    ConductivityMatproUO2
    ConstantDensityUO2
    EmissivityRelapUO2
    HeatCapacityMatproUO2
    ConstantPoissonRatioUO2
    ThermalExpansionRelapUO2
    YoungModulusMatproUO2
    DensificationFrapcon
    SwellingFrapcon
    RelocationFrapcon
    FailureUO2meltingMatpro
    PoreVelocityUO2Sens
    UO2


UPuO2
=====


Zircalloy
=========

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template: myclassinherit.rst

    DensityIAEAZy
    HeatCapacityIAEAZy
    ConductivityRelapZy
    ConstantEmissivityZy
    YoungModulusMatproZy
    ConstantPoissonRatioZy
    ThermalExpansionMatproZy
    SwellingGrowthBISONZy
    PhaseTransitionZyDynamic
    Zircalloy


HastelloyN
==========

Inconel600
==========

Steel 1515Ti
============

Molybdenum
==========

Buffer
======
