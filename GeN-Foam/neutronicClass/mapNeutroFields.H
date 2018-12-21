
tmp<volScalarField> tRhoRhok = rho*rhok;

neutroToFluid.mapTgtToSrc( TavFuel, plusEqOp<scalar>(), TfuelNeutro);
neutroToFluid.mapTgtToSrc( TavClad, plusEqOp<scalar>(), TcladNeutro);
neutroToFluid.mapTgtToSrc( tRhoRhok(), plusEqOp<scalar>(), rhoCoolNeutro);
neutroToFluid.mapTgtToSrc( thermo.T(), plusEqOp<scalar>(), TCoolNeutro);

thermo.T().correctBoundaryConditions();

tRhoRhok.clear();

if(liquidFuel)
{   
    U.correctBoundaryConditions();
    neutroToFluid.mapTgtToSrc( U, plusEqOp<vector>(), UNeutro.primitiveFieldRef());
    UNeutro.correctBoundaryConditions();

    neutroToFluid.mapTgtToSrc( porosity, plusEqOp<scalar>(), porosityNeutro.primitiveFieldRef());
    porosityNeutro.correctBoundaryConditions();

    phiNeutro = linearInterpolate(UNeutro) & neutroMesh.Sf();
}


