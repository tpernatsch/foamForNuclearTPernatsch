Thermal solver
==============

The ``thermalSubSolver`` class handles heat transfer modelling in the thermo-mechanics
module (OFFBEAT / GeN-Foam solid region).

The primary field is the temperature ``T`` (in K). A ``T`` field file must be present
at the start of the simulation to define the internal field and boundary conditions.

-------------------------------------------------------------------------------

.. rubric:: Selection in ``solverDict``

The thermal solver is selected in ``solverDict`` via the ``thermalSolver`` dictionary
(located in ``constant/`` for OFFBEAT, or in ``constant/<regionName>/`` for GeN-Foam):

.. code-block:: cpp

   thermalSolver
   {
       type solidConduction;
       // additional options (if any)
   }

-------------------------------------------------------------------------------

The thermo-mechanics module currently provides the following thermal solvers:

.. toctree::
   :maxdepth: 1
   :caption: Thermal solver models
   :glob:

   ../../../../../cppapi/generated/offbeatLib/physicsSubSolvers/thermalSubSolver/*