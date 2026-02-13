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

   foamForNuclear.porous_medium.drag.DragModel
   foamForNuclear.porous_medium.drag.DragByRegime
   foamForNuclear.porous_medium.drag.BaxiDalleDonne
   foamForNuclear.porous_medium.drag.Churchill
   foamForNuclear.porous_medium.drag.Colebrook
   foamForNuclear.porous_medium.drag.Engel
   foamForNuclear.porous_medium.drag.ModifiedEngel
   foamForNuclear.porous_medium.drag.NoKazimiFluidStructureDragModel
   foamForNuclear.porous_medium.drag.Rehme
   foamForNuclear.porous_medium.drag.ReynoldsPower
   foamForNuclear.porous_medium.drag.Autruffe
   foamForNuclear.porous_medium.drag.Bestion
   foamForNuclear.porous_medium.drag.BestionTRACE
   foamForNuclear.porous_medium.drag.NoKazimiFluidFluidDragModel
   foamForNuclear.porous_medium.drag.SchillerNaumann
   foamForNuclear.porous_medium.drag.Wallis


Heat-transfer models
====================

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   foamForNuclear.porous_medium.heat_transfer.HeatTransferModel
   foamForNuclear.porous_medium.heat_transfer.Gorenflo
   foamForNuclear.porous_medium.heat_transfer.NusseltReynoldsPrandtlPower
   foamForNuclear.porous_medium.heat_transfer.NusseltAndWall
   foamForNuclear.porous_medium.heat_transfer.NusseltWallAndHfromFMU
   foamForNuclear.porous_medium.heat_transfer.Shah
   foamForNuclear.porous_medium.heat_transfer.HeatTransferByRegime
   foamForNuclear.porous_medium.heat_transfer.NoKazimiHeatTransferModel
   foamForNuclear.porous_medium.heat_transfer.NusseltReynoldsPrandtlPowerFluidFluid
   foamForNuclear.porous_medium.heat_transfer.ConstantHeatTransfer
   foamForNuclear.porous_medium.heat_transfer.FlowEnhancementFactor
   foamForNuclear.porous_medium.heat_transfer.CobraTfFlowEnhancementFactor
   foamForNuclear.porous_medium.heat_transfer.RezkallahSimsFlowEnhancementFactor
   foamForNuclear.porous_medium.heat_transfer.SuppressionFactor
   foamForNuclear.porous_medium.heat_transfer.CobraTfSuppressionFactor
   foamForNuclear.porous_medium.heat_transfer.ChenSuppressionFactor
   foamForNuclear.porous_medium.heat_transfer.SuperpositionNucleateBoiling
   foamForNuclear.porous_medium.heat_transfer.NucleateBoilingOnsetModel
   foamForNuclear.porous_medium.heat_transfer.BasuNucleateBoilingOnsetModel
   foamForNuclear.porous_medium.heat_transfer.SubCooledBoilingFractionModel
   foamForNuclear.porous_medium.heat_transfer.SahaZuberSubCooledBoilingFractionModel
   foamForNuclear.porous_medium.heat_transfer.CriticalHeatFluxModel
   foamForNuclear.porous_medium.heat_transfer.ConstantCHFCriticalHeatFluxModel
   foamForNuclear.porous_medium.heat_transfer.LeidenfrostModel
   foamForNuclear.porous_medium.heat_transfer.GroeneveldStewartLeidenfrostModel
   foamForNuclear.porous_medium.heat_transfer.AnnularFlowModel
   foamForNuclear.porous_medium.heat_transfer.CachardLiquidAnnularFlowModel
   foamForNuclear.porous_medium.heat_transfer.CachardVapourAnnularFlowModel
   foamForNuclear.porous_medium.heat_transfer.MultiRegimeBoilingTRACECHF
   foamForNuclear.porous_medium.heat_transfer.MultiRegimeBoilingVapourTRACE


Power models
============

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   foamForNuclear.porous_medium.power_models.PowerModel
   foamForNuclear.porous_medium.power_models.FixedPower
   foamForNuclear.porous_medium.power_models.FixedTemperature
   foamForNuclear.porous_medium.power_models.FixedTemperatureFMU
   foamForNuclear.porous_medium.power_models.HeatedPin
   foamForNuclear.porous_medium.power_models.NuclearFuelPin
   foamForNuclear.porous_medium.power_models.NuclearSteadyStatePebble
   foamForNuclear.porous_medium.power_models.LumpedNuclearStructure
   foamForNuclear.porous_medium.power_models.XYPosLattice
   foamForNuclear.porous_medium.power_models.NuclearFuelFMU


Pump models
===========

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   foamForNuclear.porous_medium.Pump


Heat exchanger models
=====================

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   foamForNuclear.porous_medium.HeatExchangerModel


Regime map models
=================

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   foamForNuclear.porous_medium.regime_map.RegimeMapModel
   foamForNuclear.porous_medium.regime_map.OneParameter


Power-off criterion models
==========================

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   foamForNuclear.porous_medium.powerOffCriterionModels.PowerOffCriterionModel
   foamForNuclear.porous_medium.powerOffCriterionModels.TimerPowerOffCriterionModel
   foamForNuclear.porous_medium.powerOffCriterionModels.FieldValuePowerOffCriterionModel


Pair geometry models
====================

Mainly used in ``twoPhase`` solver.

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   foamForNuclear.porous_medium.pair_geometry.dispersion.DispersionModel
   foamForNuclear.porous_medium.pair_geometry.dispersion.ByRegime
   foamForNuclear.porous_medium.pair_geometry.dispersion.Constant
   foamForNuclear.porous_medium.pair_geometry.interfacial_area_density.InterfacialAreaDensityModel
   foamForNuclear.porous_medium.pair_geometry.interfacial_area_density.Spherical
   foamForNuclear.porous_medium.pair_geometry.contact_partition.ContactPartitionModel
   foamForNuclear.porous_medium.pair_geometry.contact_partition.ByRegime
   foamForNuclear.porous_medium.pair_geometry.contact_partition.Constant
   foamForNuclear.porous_medium.pair_geometry.PairGeometryModel


Phase change models
===================

Mainly used in ``twoPhase`` solver.

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   foamForNuclear.porous_medium.phase_change.latent_heat.LatentHeatModel
   foamForNuclear.porous_medium.phase_change.latent_heat.FinkLeibowitz
   foamForNuclear.porous_medium.phase_change.latent_heat.FromThermophysicalProperties
   foamForNuclear.porous_medium.phase_change.latent_heat.Water
   foamForNuclear.porous_medium.phase_change.saturation.SaturationModel
   foamForNuclear.porous_medium.phase_change.saturation.BrowningPotter
   foamForNuclear.porous_medium.phase_change.saturation.ConstantTemperature
   foamForNuclear.porous_medium.phase_change.saturation.Water
   foamForNuclear.porous_medium.phase_change.saturation.WaterTRACE
   foamForNuclear.porous_medium.phase_change.PhaseChangeModel
   foamForNuclear.porous_medium.phase_change.ForcedConstant
   foamForNuclear.porous_medium.phase_change.HeatDriven


Two phase drag multiplier models
================================

Mainly used in ``twoPhase`` solver.

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   foamForNuclear.porous_medium.two_phase_drag_multiplier.TwoPhaseDragMultiplierModel
   foamForNuclear.porous_medium.two_phase_drag_multiplier.LottesFlinn
   foamForNuclear.porous_medium.two_phase_drag_multiplier.LockhartMartinelli
   foamForNuclear.porous_medium.two_phase_drag_multiplier.Kaiser88
