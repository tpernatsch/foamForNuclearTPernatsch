.. _userguide_neutronics_subsolvers:

Sub-solvers
-----------

Neutronics calculations are performed by classes derived from
:ref:`neutronics.H <neutronicsModel>` that contain specific
sub-solvers. For the user, the derived classes translate into runtime selectable models. The
specific sub-solver to be used in a simulation is normally selected at the level of the application that makes use of the modules. For example, in GeN-Foam, the module is chosen in the ``regionsDict`` (see :ref:`Achieving coupled solutions <couplingGF>`).

