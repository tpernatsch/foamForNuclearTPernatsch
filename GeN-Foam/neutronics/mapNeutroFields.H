neutroToFluid.mapTgtToSrc( TavFuel, plusEqOp<scalar>(), TfuelNeutro.primitiveFieldRef());
neutroToFluid.mapTgtToSrc( TavClad, plusEqOp<scalar>(), TcladNeutro.primitiveFieldRef());
neutroToFluid.mapTgtToSrc( thermo.rho(), plusEqOp<scalar>(), rhoCoolNeutro.primitiveFieldRef());
neutroToFluid.mapTgtToSrc( thermo.T(), plusEqOp<scalar>(), TCoolNeutro.primitiveFieldRef());
if(liquidFuel)
{	

	U.correctBoundaryConditions();
	neutroToFluid.mapTgtToSrc( U, plusEqOp<vector>(), UNeutro.primitiveFieldRef());
	UNeutro.correctBoundaryConditions();

	neutroToFluid.mapTgtToSrc( porosity, plusEqOp<scalar>(), porosityNeutro.primitiveFieldRef());
	porosityNeutro.correctBoundaryConditions();

	phiNeutro = linearInterpolate(UNeutro) & neutroMesh.Sf();
}


