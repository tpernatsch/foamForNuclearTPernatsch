Mechanics solver
================

The ``mechanicsSubSolver`` class handles the solution of the solid mechanics problem
in the thermo-mechanics module (OFFBEAT / GeN-Foam solid region), i.e. the balance of
linear momentum.

The primary field is the displacement. Depending on the solver formulation this is:

* ``D``  : total displacement (m)
* ``DD`` : incremental displacement (m)

The corresponding displacement field file must be present at the start of the simulation
to define the internal field and boundary conditions.

-------------------------------------------------------------------------------

.. rubric:: Selection in ``solverDict``

The mechanics solver is selected in ``solverDict`` via the ``mechanicsSolver`` dictionary
(located in ``constant/`` for OFFBEAT, or in ``constant/<regionName>/`` for GeN-Foam):

.. code-block:: cpp

   mechanicsSolver
   {
       type smallStrain;
       // additional options (if any)
   }

.. note::

   The mechanics solver often benefits from the use of the
   :doc:`multi-material interface correction <multi_material_correction>`
   for improved accuracy across material boundaries.

-------------------------------------------------------------------------------

.. warning::

   When using an incremental or updated formulation, boundary and initial conditions
   **must** be defined for ``DD`` (not for ``D``). In such cases, ``D`` does not need to
   be present in the initial time folder, as it is reconstructed during the simulation.
   The linear solver settings (and residual controls) must then also be provided for
   ``DD`` in the system solution dictionaries.

-------------------------------------------------------------------------------

The thermo-mechanics module currently provides the following mechanics solvers:

.. toctree::
   :maxdepth: 1
   :caption: Mechanics solver models
   :glob:

   ../../../../../cppapi/generated/offbeatLib/physicsSubSolvers/mechanicsSubSolver/*