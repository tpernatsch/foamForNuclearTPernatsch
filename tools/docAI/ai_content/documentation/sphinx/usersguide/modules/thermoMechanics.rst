.. _userguide_thermomechanics:

Thermal-mechanics
=================

GeN-Foam can currently use two different thermomechanics solvers. One is a
simple linear elasticity solver, and the other is an extended thermomechanics
solver derived from OFFBEAT [:ref:`SCOLARO2020110416 <SCOLARO2020110416>`]. The former
is presented here, whereas the latter is described in detail in the
`OFFBEAT Documentation <https://foam-for-nuclear.gitlab.io/offbeat/>`_.

.. note ::

    In OFFBEAT, what is named *solverDict* is called *thermomechanicalProperties* in GeN-Foam, to ensure continuity between the two mechanics solvers.


Various properties
------------------

The *thermoMechanicalProperties* dictionary
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The *thermoMechanicalProperties* dictionary can be found under
*constant/thermoMechanicalRegion*. It allows you to define the
thermo-mechanical properties of structures, subdivided according to the
cellZones of the *thermoMechanicalRegion* mesh.

A detailed, commented example is available in the tutorial
`3D_SmallESFR <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/3D_SmallESFR/extendedThermoMechanics/constant/thermoMechanicalRegion/thermoMechanicalProperties>`_.


Initial and boundary conditions
-------------------------------

In addition to the standard ones available in OpenFOAM, GeN-Foam includes a
*tractionDisplacement* boundary condition. This allows you to prescribe a
pressure or a traction on a boundary. Work is ongoing to adapt to GeN-Foam
the contact boundary condition from the OFFBEAT fuel behavior solver
[:ref:`SCOLARO2020110416 <SCOLARO2020110416>`].


Discretization and solution
---------------------------

Details for the discretization and solution of the equations are handled in a
standard OpenFOAM way, through the *fvSolution* and *fvSchemes* dictionaries
in *constant/thermoMechanicalRegion*.


Mesh deformation
----------------

The new structure of GeN-Foam allows you to deform each mesh based on any
user-defined vectorial field. In particular, it may be of interest to solve
the neutron transport on a deformed mesh. To do this, the user needs to
ensure that the *meshDisp* field is mapped from the mechanics to the
neutronics.

Moreover, in *constant/regionsDict*, you can specify for each physics
whether the mesh needs to be deformed, and which field is used to compute
the deformation. An example is provided in the
`3D_SmallESFR <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/3D_SmallESFR/extendedThermoMechanics/system/controlDict>`_
tutorial.