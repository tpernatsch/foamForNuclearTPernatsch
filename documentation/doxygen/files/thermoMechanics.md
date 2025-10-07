
# Thermal-mechanics {#TM}

GeN-Foam can currently use two different thermomechanics solvers. One is a simple linear elasticity solver and the other is an extended thermomechanics solver derived from OFFBEAT. The former will be presented here, whereas the latter is described in details in the OFFBEAT guide.

N.B. what in OFFBEAT is named *solverDict*, in GeN-Foam has been renamed *thermomechanicalProperties*, for continuity between the two mechanics solvers.

## Various properties

<div class="border-box">
<b>The *thermoMechanicalProperties* dictionary</b>

The *thermoMechanicalProperties* dictionary can be found under *constant/thermoMechanicalRegion* and allow to define the thermo-mechanical properties of structures, subdivided according to the cellZones of the thermoMechanicalRegion mesh.
<br><br>One can find a detailed, commented example in the tutorial
[3D_SmallESFR](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/reactorCases/3D_SmallESFR_NewSolverVerification/newSolver/constant/thermoMechanicalRegion/thermoMechanicalProperties).
</div>
<br>


## Initial and boundary conditions

Besides the standard ones available in OpenFOAM, GeN-Foam includes a *tractionDisplacement* boundary condition that allows to set a pressure or a traction on a boundary. Work is ongoing to adapt to GeN-Foam the contact boundary condition of the OFFBEAT fuel behavior solver \cite SCOLARO2020110416.


## Discretization and solution

Details for discretization and solution of equations are handled in a standard OpenFOAM way, i.e., through the *fvSolution* and *fvSchemes* dictionaries in *constant/thermoMechanicalRegion*.


## Mesh deformation

The new structure of GeN-Foam allows to deform each mesh based on any user-defined vectorial field. In particular, it might be of interest to solve the neutron transport on a deformed mesh. To do this, the user needs to make sure to map the *meshDisp* field from the mechanics to the neutronics. Moreover, in the *constant/multiRegionCouplingDict*, one can specify for each physics if the mesh needs to be deformed and based on which field the deformation is computed. An example of this is found in the [3D_SmallESFR](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/reactorCases/3D_SmallESFR_NewSolverVerification/newSolver/system/controlDict) tutorial.
