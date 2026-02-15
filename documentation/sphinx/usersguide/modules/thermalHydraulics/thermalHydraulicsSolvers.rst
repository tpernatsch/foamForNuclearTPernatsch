.. _userguide_thermalhydraulics_subsolvers:

Sub-solvers
-----------

Thermal-hydraulics calculations are performed by classes derived from
:ref:`thermalHydraulicsModel.H <thermalHydraulicsModel>` that contain specific
sub-solvers. For the user, the derived classes translate into runtime selectable models. The
specific sub-solver to be used in a simulation is normally selected at the level of the application that makes use of the modules. For example, in GeN-Foam, the module is chosen in the ``regionsDict`` (see :ref:`Achieving coupled solutions <couplingGF>`). Available sub-solvers include:

.. toctree::
   :maxdepth: 1
   :glob:

    ../../../sphinx/cppapi/generated/modules/thermalHydraulics/**

.. raw:: html

   <br><br>


.. See refs:
    :ref:`RADMAN2021111178 <RADMAN2021111178>`
    :ref:`RADMAN2021111422 <RADMAN2021111422>`



