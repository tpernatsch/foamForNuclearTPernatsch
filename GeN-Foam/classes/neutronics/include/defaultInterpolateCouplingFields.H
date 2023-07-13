neutroToFluid.mapTgtToSrc(*TFuelOrig_, plusEqOp<scalar>(), TFuel_);
neutroToFluid.mapTgtToSrc(*TCladOrig_, plusEqOp<scalar>(), TClad_);
neutroToFluid.mapTgtToSrc(*TCoolOrig_, plusEqOp<scalar>(), TCool_);
neutroToFluid.mapTgtToSrc(*rhoCoolOrig_, plusEqOp<scalar>(), rhoCool_);
if (liquidFuel_)
{
    neutroToFluid.mapTgtToSrc(*UOrig_, plusEqOp<vector>(), UPtr_());
    neutroToFluid.mapTgtToSrc
    (
        *alphaOrig_, 
        plusEqOp<scalar>(), 
        alphaPtr_()
    );
    neutroToFluid.mapTgtToSrc
    (
        *alphatOrig_, 
        plusEqOp<scalar>(), 
        alphatPtr_()
    );
    neutroToFluid.mapTgtToSrc(*muOrig_, plusEqOp<scalar>(), muPtr_());
    phiPtr_() = fvc::flux(UPtr_());
    volScalarField diffCoeffOrig
    (
        (
            *alphatOrig_ 
        +   *muOrig_/xs_.ScNo()
        )/(*rhoCoolOrig_)
    ); 
    neutroToFluid.mapTgtToSrc
    (
        diffCoeffOrig, 
        plusEqOp<scalar>(), 
        diffCoeffPrecPtr_()
    );

    UPtr_().correctBoundaryConditions();
    alphaPtr_().correctBoundaryConditions();
    alphatPtr_().correctBoundaryConditions();
    diffCoeffPrecPtr_().correctBoundaryConditions();
}

//- Interpolate T from TM solver
// This can be a bit confusing: we need to paramtrize XS based on this
// temperature, but XS can only be paramterized based on 
// TFuel_, TClad_, rhoCool_, TCool_, disp_
// I am using TFuel and this os purely arbitrary. It could have been any
// of the other temperatures. However, this means that one will have to
// provide a nuclearDataFuelTemp file to parametrize the XS from a solid
// structure
neutroToMech.mapTgtToSrc(*TFuelOrigMech_, plusEqOp<scalar>(), TFuel_);

TFuel_.correctBoundaryConditions();
TClad_.correctBoundaryConditions();
TCool_.correctBoundaryConditions();
rhoCool_.correctBoundaryConditions();


