.. _pythonapi_base:

--------------------------------------------
:mod:`foamForNuclear` -- Basic Functionality
--------------------------------------------

.. contents:: Table of Contents
    :local:


Pre-processing
==============

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   foamForNuclear.mesh.dicts.CreateBafflesDict
   foamForNuclear.mesh.dicts.DecomposeParDict
   foamForNuclear.preprocessing.SetFieldsDict
   foamForNuclear.mesh.dicts.TopoSetDict

More information on meshing in :mod:`foamForNuclear.mesh`.


Simulation Settings
===================

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   foamForNuclear.control.ControlDict
   foamForNuclear.solvers.Solver
   foamForNuclear.numerics.fvSchemes
   foamForNuclear.numerics.fvSolution
   foamForNuclear.quadratureSet.QuadratureSet
   foamForNuclear.transport.TransportProperties
   foamForNuclear.turbulence.TurbulenceProperties


Fields Specifications
=====================

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   foamForNuclear.timeFolder.TimeFolder
   foamForNuclear.fields.Field
   foamForNuclear.fields.Dimension
   foamForNuclear.fields.ReducedDimension


Solvers
=======

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   foamForNuclear.solvers.Solver
   foamForNuclear.solvers.ThermalHydraulicsSolver
   foamForNuclear.solvers.CompressibleInterFoam
   foamForNuclear.solvers.NeutronicsSolver
   foamForNuclear.solvers.OffbeatSolver

More information on thermal-hydraulics sub-models in
:doc:`thermalHydraulicsModels`.


Nuclear Data and Neutronics dictionaries
========================================

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   foamForNuclear.nuclearData.NuclearDataZone
   foamForNuclear.nuclearData.NuclearDataState
   foamForNuclear.nuclearData.NuclearData
   foamForNuclear.nuclearData.PointKineticsData
   foamForNuclear.quadratureSet.QuadratureSet
   foamForNuclear.ExternalSource
   foamForNuclear.nuclearData.ControlRodMove
   foamForNuclear.nuclearData.ControlRodMovement


Thermal-mechanics and Fuel Performance
======================================

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   foamForNuclear.offbeat_lib.ThermoMechanicsCouplingOptions
   foamForNuclear.offbeat_lib.GlobalOptions
   foamForNuclear.offbeat_lib.ThermalSolverOptions
   foamForNuclear.offbeat_lib.ReadTemperatureThermalSolverOptions
   foamForNuclear.offbeat_lib.SolidConductionThermalSolverOptions
   foamForNuclear.offbeat_lib.MechanicsSolverOptions
   foamForNuclear.offbeat_lib.SmallStrainMechanicsSolverOptions
   foamForNuclear.offbeat_lib.SmallStrainIncrementalUpdatedMechanicsSolverOptions
   foamForNuclear.offbeat_lib.LargeStrainTotLagMechanicsSolverOptions
   foamForNuclear.offbeat_lib.LargeStrainUpdLagMechanicsSolverOptions
   foamForNuclear.offbeat_lib.NeutronicsSolverOptions
   foamForNuclear.offbeat_lib.DiffusionNeutronicsSolverOptions
   foamForNuclear.offbeat_lib.ElementTransportSolverOptions
   foamForNuclear.offbeat_lib.ByListElementTransportSolverOptions
   foamForNuclear.offbeat_lib.burnup.BurnupOptions
   foamForNuclear.offbeat_lib.ConstantBurnupOptions
   foamForNuclear.offbeat_lib.FromPowerBurnupOptions
   foamForNuclear.offbeat_lib.LassmannBurnupOptions
   foamForNuclear.offbeat_lib.FgrOptions
   foamForNuclear.offbeat_lib.SciantixFgrOptions
   foamForNuclear.offbeat_lib.GapGasOptions
   foamForNuclear.offbeat_lib.TrisoGapGasOptions
   foamForNuclear.offbeat_lib.FrapconGapGasOptions
   foamForNuclear.offbeat_lib.RheologyOptions
   foamForNuclear.offbeat_lib.ByMaterialRheologyOptions
   foamForNuclear.offbeat_lib.StressAnalysis
   foamForNuclear.offbeat_lib.HeatSourceOptions
   foamForNuclear.offbeat_lib.TimeDependentLhgrHeatSourceOptions
   foamForNuclear.offbeat_lib.TimeDependentVhgrHeatSourceOptions
   foamForNuclear.offbeat_lib.FastFluxOptions
   foamForNuclear.offbeat_lib.TimeDependentFastFluxOptions
   foamForNuclear.offbeat_lib.CorrosionOptions
   foamForNuclear.offbeat_lib.SliceMapperOptions
   foamForNuclear.offbeat_lib.ByMaterialSliceMapperOptions
   foamForNuclear.offbeat_lib.AutoAxialSliceMapperOptions
   foamForNuclear.offbeat_lib.ByPelletSliceMapperOptions


Coupling
========

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   foamForNuclear.coupling.FieldTransfer
   foamForNuclear.coupling.MultiPhysicsLoop
   foamForNuclear.coupling.CHTLoop
   foamForNuclear.coupling.FSILoop
   foamForNuclear.coupling.Coupling
   foamForNuclear.coupling.ExternalCouplingDict


Post-processing
===============

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   foamForNuclear.functions.FunctionObject
   foamForNuclear.functions.SurfaceFieldValue
   foamForNuclear.functions.VolFieldValue
   foamForNuclear.functions.FieldMinMax
   foamForNuclear.functions.MapFields
   foamForNuclear.functions.Probes
   foamForNuclear.functions.MassFlow
   foamForNuclear.functions.TBulk
   foamForNuclear.functions.FMUSimulator
   foamForNuclear.functions.FunctionObjects


Mesh Motion
===========

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   foamForNuclear.mesh.dicts.DynamicMeshDict
   foamForNuclear.mesh.dicts.MotionDict
   foamForNuclear.mesh.dicts.LinearMotion
   foamForNuclear.mesh.dicts.OscillatingLinearMotion
   foamForNuclear.mesh.dicts.RotatingMotion
   foamForNuclear.mesh.dicts.OscillatingRotatingMotion
   foamForNuclear.mesh.dicts.DiffusivityMotion
   foamForNuclear.mesh.dicts.UniformDiffusivityMotion
   foamForNuclear.mesh.dicts.DirectionalDiffusivityMotion
   foamForNuclear.mesh.dicts.MotionDirectionalDiffusivityMotion
   foamForNuclear.mesh.dicts.QuadraticDiffusivityMotion
   foamForNuclear.mesh.dicts.FileDiffusivityMotion

More information on mesh motion in :ref:`Mesh motion <pythonapi_mesh_motion>`.


Common objects
==============

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   foamForNuclear.timeProfile.TimeProfile
