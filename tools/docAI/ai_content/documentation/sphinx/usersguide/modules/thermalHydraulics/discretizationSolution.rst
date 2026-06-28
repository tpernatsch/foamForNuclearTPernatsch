.. _userguide_thermalhydraulics_discretizationSolution:


Discretization and solution
----------------------------

Details for the discretization and solution of the equations are handled in a
standard OpenFOAM way, i.e., through the *fvSolution* and *fvSchemes*
dictionaries in *system/[nameOfTheFluidRegion]*. Note that certain additions
may have been included depending on the specific
(:ref:`solver <userguide_thermalhydraulics_subsolvers>`).