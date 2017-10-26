
tmp<volScalarField> tRhoRhok = thermo.rho()*rhok;

neutroToFluid.mapTgtToSrc( TavFuel, plusEqOp<scalar>(), TfuelNeutro.primitiveFieldRef());
neutroToFluid.mapTgtToSrc( TavClad, plusEqOp<scalar>(), TcladNeutro.primitiveFieldRef());
neutroToFluid.mapTgtToSrc( tRhoRhok(), plusEqOp<scalar>(), rhoCoolNeutro.primitiveFieldRef());
neutroToFluid.mapTgtToSrc( thermo.T(), plusEqOp<scalar>(), TCoolNeutro.primitiveFieldRef());

tRhoRhok.clear();

TfuelNeutro.correctBoundaryConditions();
TcladNeutro.correctBoundaryConditions();
rhoCoolNeutro.correctBoundaryConditions();
TCoolNeutro.correctBoundaryConditions();

if(liquidFuel)
{   
    U.correctBoundaryConditions();
    neutroToFluid.mapTgtToSrc( U, plusEqOp<vector>(), UNeutro.primitiveFieldRef());
    UNeutro.correctBoundaryConditions();

    neutroToFluid.mapTgtToSrc( porosity, plusEqOp<scalar>(), porosityNeutro.primitiveFieldRef());
    porosityNeutro.correctBoundaryConditions();

    phiNeutro = linearInterpolate(UNeutro) & neutroMesh.Sf();
}


