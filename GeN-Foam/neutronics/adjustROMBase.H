
if(adjustROM)
{
    Info << "Adjusting POD base" << endl ;

    Info << "ROMAdjustCounter = " << ROMAdjustCounter << endl ;

    if( ROMAdjustCounter == 1 )//first adjustement in step - base clean up
    {

            base_fluxStar[energyI].resize(1);
            base_fluxStar[energyI].set
            (
                0,
                new volScalarField 
                (
                    IOobject
                    (
                        "_base_fluxStar"+Foam::name(energyI)+"_"+Foam::name(0),
                        runTime.timeName(),
                        neutroMesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    fluxStar[energyI]
                )
            );
            base_fluxStar[energyI][0].primitiveFieldRef() = fluxStar[energyI].internalField() / Foam::pow(Foam::POD::projection(fluxStar[energyI],fluxStar[energyI]),0.5);
            base_fluxStar[energyI][0].correctBoundaryConditions();
    }
    else
    {

        int initBaseSize = base_fluxStar[energyI].size();
        base_fluxStar[energyI].resize(initBaseSize+1);
        base_fluxStar[energyI].set
        (
            initBaseSize,
            new volScalarField 
            (
                IOobject
                (
                    "_base_fluxStar"+Foam::name(energyI)+"_"+Foam::name(initBaseSize),
                    runTime.timeName(),
                    neutroMesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fluxStar[energyI]
            )
        );
        for(int baseNoI = 0; baseNoI < initBaseSize ; baseNoI++ )
        {
            base_fluxStar[energyI][initBaseSize]
             -= Foam::POD::projection(fluxStar[energyI],base_fluxStar[energyI][baseNoI])*base_fluxStar[energyI][baseNoI];
        }
        base_fluxStar[energyI][initBaseSize].correctBoundaryConditions();
        base_fluxStar[energyI][initBaseSize] /= Foam::pow(Foam::POD::projection(base_fluxStar[energyI][initBaseSize],base_fluxStar[energyI][initBaseSize]),0.5);
        base_fluxStar[energyI][initBaseSize].correctBoundaryConditions();   

        Info << "PODBaseSize energy " << energyI << ": " << base_fluxStar[energyI].size() << endl;
        Info << "Done" << endl;

        #include "buildMatrixPOD.H"

    }

}
else
{
    Info << "Not adjusting POD base" << endl ;
}



