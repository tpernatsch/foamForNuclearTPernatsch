Neutronics solver
=================

The ``neutronicsSubSolver`` class handles the computation of neutron flux
within fuel pellets in the thermo-mechanics module
(OFFBEAT / GeN-Foam solid region).

While fuel performance codes typically do not solve neutronics equations,
such calculations can become necessary when radial or symmetric flux
profiles cannot be assumed, as is commonly done in classical 1.5D
simulations. In those cases, solving a neutron transport approximation
directly in the fuel domain allows spatial variations of the neutron
flux to be captured.

In future developments, the neutron flux could also serve as input for
depletion calculations if a full depletion module is introduced.

-------------------------------------------------------------------------------

.. rubric:: Selection in ``solverDict``

The neutronics solver is selected in ``solverDict`` via the
``neutronicsSolver`` dictionary (located in ``constant/`` for OFFBEAT,
or in ``constant/<regionName>/`` for GeN-Foam):

.. code-block:: cpp

   neutronicsSolver
   {
       type diffusion;
       // additional options (if any)
   }

If the ``neutronicsSolver`` dictionary is not present in ``solverDict``,
the neutronics solution is deactivated and no additional fields are created.

-------------------------------------------------------------------------------

The thermo-mechanics module currently provides the following neutronics
solvers:

.. toctree::
   :maxdepth: 1
   :caption: Neutronics solver models
   :glob:

   ../../../../../cppapi/generated/offbeatLib/physicsSubSolvers/neutronicsSubSolver/*