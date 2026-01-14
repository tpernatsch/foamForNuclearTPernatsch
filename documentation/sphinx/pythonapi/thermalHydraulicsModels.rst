.. _pythonapi_thermalHydraulicsModels:

-------------------------
Thermal-hydraulics models
-------------------------

.. contents:: Table of Contents
    :local:

Drag models
===========

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   foamForNuclear.DragModel
   foamForNuclear.DragByRegime
   foamForNuclear.BaxiDalleDonne
   foamForNuclear.Churchill
   foamForNuclear.Colebrook
   foamForNuclear.Engel
   foamForNuclear.ModifiedEngel
   foamForNuclear.NoKazimiFluidStructureDragModel
   foamForNuclear.Rehme
   foamForNuclear.ReynoldsPower
   foamForNuclear.Autruffe
   foamForNuclear.Bestion
   foamForNuclear.BestionTRACE
   foamForNuclear.NoKazimiFluidFluidDragModel
   foamForNuclear.SchillerNaumann
   foamForNuclear.Wallis


Heat-transfer models
====================

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   foamForNuclear.HeatTransferModel
   foamForNuclear.Gorenflo
   foamForNuclear.NusseltReynoldsPrandtlPower
   foamForNuclear.NusseltAndWall
   foamForNuclear.NusseltWallAndHfromFMU
   foamForNuclear.Shah
   foamForNuclear.HeatTransferByRegime
   foamForNuclear.NoKazimiHeatTransferModel
   foamForNuclear.NusseltReynoldsPrandtlPowerFluidFluid
   foamForNuclear.ConstantHeatTransfer
   foamForNuclear.FlowEnhancementFactor
   foamForNuclear.CobraTfFlowEnhancementFactor
   foamForNuclear.RezkallahSimsFlowEnhancementFactor
   foamForNuclear.SuppressionFactor
   foamForNuclear.CobraTfSuppressionFactor
   foamForNuclear.ChenSuppressionFactor
   foamForNuclear.SuperpositionNucleateBoiling
   foamForNuclear.NucleateBoilingOnsetModel
   foamForNuclear.BasuNucleateBoilingOnsetModel
   foamForNuclear.SubCooledBoilingFractionModel
   foamForNuclear.SahaZuberSubCooledBoilingFractionModel
   foamForNuclear.CriticalHeatFluxModel
   foamForNuclear.ConstantCHFCriticalHeatFluxModel
   foamForNuclear.LeidenfrostModel
   foamForNuclear.GroeneveldStewartLeidenfrostModel
   foamForNuclear.AnnularFlowModel
   foamForNuclear.CachardLiquidAnnularFlowModel
   foamForNuclear.CachardVapourAnnularFlowModel
   foamForNuclear.MultiRegimeBoilingTRACECHF
   foamForNuclear.MultiRegimeBoilingVapourTRACE


Power models
============

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   foamForNuclear.PowerModel
   foamForNuclear.FixedPower
   foamForNuclear.FixedTemperature
   foamForNuclear.FixedTemperatureFMU
   foamForNuclear.HeatedPin
   foamForNuclear.NuclearFuelPin
   foamForNuclear.NuclearSteadyStatePebble
   foamForNuclear.LumpedNuclearStructure
   foamForNuclear.XYPosLattice
   foamForNuclear.NuclearFuelFMU


Pump models
===========

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   foamForNuclear.Pump


Heat exchanger models
=====================

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   foamForNuclear.HeatExchangerModel


Regime map models
=================

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   foamForNuclear.RegimeMapModel
   foamForNuclear.RegimeMapOneParameter


Power-off criterion models
==========================

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   foamForNuclear.PowerOffCriterionModel
   foamForNuclear.TimerPowerOffCriterionModel
   foamForNuclear.FieldValuePowerOffCriterionModel


Pair geometry models
====================

Mainly used in ``twoPhase`` solver.

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   foamForNuclear.DispersionModel
   foamForNuclear.ByRegimeDispersionModel
   foamForNuclear.ConstantDispersionModel
   foamForNuclear.InterfacialAreaDensityModel
   foamForNuclear.SphericalInterfacialAreaDensityModel
   foamForNuclear.ContactPartitionModel
   foamForNuclear.ByRegimeContactPartitionModel
   foamForNuclear.ConstantContactPartitionModel
   foamForNuclear.PairGeometryModel


Phase change models
===================

Mainly used in ``twoPhase`` solver.

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   foamForNuclear.LatentHeatModel
   foamForNuclear.FinkLeibowitzLatentHeat
   foamForNuclear.FromThermophysicalPropertiesLatentHeat
   foamForNuclear.WaterLatentHeat
   foamForNuclear.SaturationModel
   foamForNuclear.BrowningPotterSaturationModel
   foamForNuclear.ConstantTemperature
   foamForNuclear.WaterSaturationModel
   foamForNuclear.WaterTRACESaturationModel
   foamForNuclear.PhaseChangeModel
   foamForNuclear.ForcedConstant
   foamForNuclear.HeatDriven


Two phase drag multiplier models
================================

Mainly used in ``twoPhase`` solver.

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   foamForNuclear.TwoPhaseDragMultiplierModel
   foamForNuclear.LottesFlinn
   foamForNuclear.LockhartMartinelli
   foamForNuclear.Kaiser88
