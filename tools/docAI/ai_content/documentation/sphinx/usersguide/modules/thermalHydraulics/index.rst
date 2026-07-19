.. _userguide_thermalhydraulics:


Thermal-hydraulics
==================

This section introduces how foamForNuclear sets up and solves thermal-hydraulics problems in OpenFOAM-based workflows, from model specification to the final discretized equations. It begins with the thermal-hydraulics sub-solvers, then explains how initial and boundary conditions define the physical state and constraints before any numerical solution is attempted. The subsequent subsections cover the discretization and solution strategy, followed by the key configuration dictionaries—`g`, `thermophysicalProperties`, `turbulenceProperties`, and `phaseProperties`—which supply gravity, material and transport properties, turbulence modeling inputs, and (when applicable) phase-related parameters that the sub-solvers use during assembly and iteration.

.. toctree::
   :numbered:
   :maxdepth: 1

   thermalHydraulicsSolvers
   initialAndBC
   discretizationSolution
   g
   thermophysicalProperties
   turbulenceProperties
   phaseProperties/index

