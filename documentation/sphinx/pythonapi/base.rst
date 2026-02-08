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

   foamForNuclear.CreateBafflesDict
   foamForNuclear.DecomposeParDict
   foamForNuclear.SetFieldsDict
   foamForNuclear.TopoSetDict

More information on meshing in :mod:`foamForNuclear.mesh`.


Simulation Settings
===================

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   foamForNuclear.ControlDict
   foamForNuclear.Solver
   foamForNuclear.fvSchemes
   foamForNuclear.fvSolution
   foamForNuclear.PhaseProperties
   foamForNuclear.QuadratureSet
   foamForNuclear.TransportProperties
   foamForNuclear.TurbulenceProperties


Fields Specifications
=====================

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   foamForNuclear.TimeFolder
   foamForNuclear.Field
   foamForNuclear.Dimension
   foamForNuclear.ReducedDimension


Solvers
=======

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   foamForNuclear.Solver
   foamForNuclear.ThermalHydraulicsSolver
   foamForNuclear.CompressibleInterFoam
   foamForNuclear.NeutronicsSolver
   foamForNuclear.OffbeatSolver

More information on thermal-hydraulics sub-models in
:doc:`thermalHydraulicsModels`.


Nuclear Data and Neutronics dictionaries
========================================

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   foamForNuclear.NuclearDataZone
   foamForNuclear.NuclearDataState
   foamForNuclear.NuclearData
   foamForNuclear.PointKineticsData
   foamForNuclear.QuadratureSet
   foamForNuclear.ExternalSource
   foamForNuclear.ControlRodMove


Thermal-mechanics and Fuel Performance
======================================

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   foamForNuclear.ThermoMechanicsCouplingOptions
   foamForNuclear.GlobalOptions
   foamForNuclear.ThermalSolverOptions
   foamForNuclear.ReadTemperatureThermalSolverOptions
   foamForNuclear.SolidConductionThermalSolverOptions
   foamForNuclear.MechanicsSolverOptions
   foamForNuclear.SmallStrainMechanicsSolverOptions
   foamForNuclear.SmallStrainIncrementalUpdatedMechanicsSolverOptions
   foamForNuclear.LargeStrainTotLagMechanicsSolverOptions
   foamForNuclear.LargeStrainUpdLagMechanicsSolverOptions
   foamForNuclear.NeutronicsSolverOptions
   foamForNuclear.DiffusionNeutronicsSolverOptions
   foamForNuclear.ElementTransportSolverOptions
   foamForNuclear.ByListElementTransportSolverOptions
   foamForNuclear.BurnupOptions
   foamForNuclear.ConstantBurnupOptions
   foamForNuclear.FromPowerBurnupOptions
   foamForNuclear.LassmannBurnupOptions
   foamForNuclear.FgrOptions
   foamForNuclear.SciantixFgrOptions
   foamForNuclear.GapGasOptions
   foamForNuclear.TrisoGapGasOptions
   foamForNuclear.FrapconGapGasOptions
   foamForNuclear.RheologyOptions
   foamForNuclear.ByMaterialRheologyOptions
   foamForNuclear.StressAnalysis
   foamForNuclear.HeatSourceOptions
   foamForNuclear.TimeDependentLhgrHeatSourceOptions
   foamForNuclear.TimeDependentVhgrHeatSourceOptions
   foamForNuclear.FastFluxOptions
   foamForNuclear.TimeDependentFastFluxOptions
   foamForNuclear.CorrosionOptions
   foamForNuclear.SliceMapperOptions
   foamForNuclear.ByMaterialSliceMapperOptions
   foamForNuclear.AutoAxialSliceMapperOptions
   foamForNuclear.ByPelletSliceMapperOptions


Coupling
========

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   foamForNuclear.FieldTransfer
   foamForNuclear.MultiPhysicsLoop
   foamForNuclear.CHTLoop
   foamForNuclear.FSILoop
   foamForNuclear.Coupling
   foamForNuclear.ExternalCouplingDict


Post-processing
===============

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   foamForNuclear.FunctionObject
   foamForNuclear.SurfaceFieldValue
   foamForNuclear.VolFieldValue
   foamForNuclear.FieldMinMax
   foamForNuclear.MapFields
   foamForNuclear.Probes
   foamForNuclear.MassFlow
   foamForNuclear.TBulk
   foamForNuclear.FMUSimulator
   foamForNuclear.FunctionObjects


Mesh Motion
===========

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   foamForNuclear.DynamicMeshDict
   foamForNuclear.MotionDict
   foamForNuclear.LinearMotion
   foamForNuclear.OscillatingLinearMotion
   foamForNuclear.RotatingMotion
   foamForNuclear.OscillatingRotatingMotion
   foamForNuclear.DiffusivityMotion
   foamForNuclear.UniformDiffusivityMotion
   foamForNuclear.DirectionalDiffusivityMotion
   foamForNuclear.MotionDirectionalDiffusivityMotion
   foamForNuclear.QuadraticDiffusivityMotion
   foamForNuclear.FileDiffusivityMotion

More information on mesh motion in :ref:`Mesh motion <pythonapi_mesh_motion>`.


Common objects
==============

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   foamForNuclear.TimeProfile
