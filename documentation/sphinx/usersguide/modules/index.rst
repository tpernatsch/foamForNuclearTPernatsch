.. _modules:

---------------
Physics modules
---------------

Within foamForNuclear, the **modules** are those parts of the code that solves for specific physics. Typically, modules offer multiple run-time selectable sub-solvers (e.g., single- and two-phase sub-solvers or diffusion and SP3 sub-solvers), as well as extensive sets of materials properties, correlations and behavioral models. 

   - Users can employ all modules within GeN-Foam while OFFBEAT is a specialized solver that makes an advanced use of a single thermoMechanics module. 
   - Developers can: use modules to create new applications instead of GeN-Foam and OFFBEAT; develop new modules by using the existing ones as templates; or expand existing modules with new sub-solvers, meaterials, correlations, etc. 

Form a programming perspective, modules are hierachical C++ structures where a parent class is used as a standard interface for using the module within applications. That parent class (e.g., neutronics) can be used to derive specialized classes (e.g., diffusion or SP3 ) and both partent and derived classes can make use of specialized libraries (e.g., cross sections). 



.. toctree::
   :maxdepth: 3

   thermalHydraulics/index
   neutronics/index
   thermoMechanics
   openfoamImportedSolvers

