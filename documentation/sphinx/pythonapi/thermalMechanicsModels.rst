.. _pythonapi_thermalMechanicsModels:

.. module:: foamForNuclear.offbeat_lib.materials

--------------------------------------------------------------------------
:mod:`foamForNuclear.offbeat_lib.materials` -- Thermo-mechanical Materials
--------------------------------------------------------------------------

.. contents:: Table of Contents
    :local:


Base and abstract thermo-mechanical class
=========================================

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template: myclassinherit.rst

    properties.BaseThermomechanicalMaterial
    properties.ThermomechanicalPropertyModel
    properties.conductivity.ConductivityModel
    properties.density.DensityModel
    properties.emissivity.EmissivityModel
    properties.heat_capacity.HeatCapacityModel
    properties.poisson_ratio.PoissonRatioModel
    properties.thermal_expansion.ThermalExpansionModel
    properties.young_modulus.YoungModulusModel
    behaviour.densification.DensificationModel
    behaviour.swelling.SwellingModel
    behaviour.phase_transition.PhaseTransitionModel
    behaviour.relocation.RelocationModel
    behaviour.failure.FailureModel
    properties.PoreVelocityModel
    laws.yield_stress.YieldStress
    laws.yield_stress.Hardening
    laws.yield_stress.Constant
    laws.yield_stress.FRAPTRAN
    laws.creep.Creep
    laws.creep.Limback
    RheologyConstitutiveLaw
    ElasticityRheologyModel
    MisesPlasticityRheologyModel
    MisesPlasticCreepRheologyModel
    Material
    FuelMaterial
    laws.ConstantMechanicalLaw


Constant Material
=================

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template: myclassinherit.rst

    Constant


UO2
===

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template: myclassinherit.rst

    properties.conductivity.UO2Matpro
    properties.density.UO2Constant
    properties.emissivity.UO2Relap
    properties.heat_capacity.UO2Matpro
    properties.poisson_ratio.UO2Constant
    properties.thermal_expansion.UO2Relap
    properties.young_modulus.UO2Matpro
    behaviour.densification.UO2Frapcon
    behaviour.failure.UO2Matpro
    behaviour.relocation.UO2Frapcon
    behaviour.swelling.UO2Frapcon
    laws.creep.UO2Frapcon
    UO2


UPuO2
=====


Zircaloy
========

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template: myclassinherit.rst

    properties.density.ZircaloyIaea
    properties.heat_capacity.ZircaloyMatpro
    properties.heat_capacity.ZircaloyIaea
    properties.conductivity.ZircaloyRelap
    properties.emissivity.ZircaloyConstant
    properties.young_modulus.ZircaloyMatpro
    properties.poisson_ratio.ZircaloyConstant
    properties.poisson_ratio.ZircaloyMatpro
    properties.thermal_expansion.ZircaloyMatpro
    behaviour.swelling.ZircaloyBison
    behaviour.swelling.ZircaloyMatpro
    behaviour.phase_transition.ZircaloyDynamic
    Zircaloy


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
