// Nuclear data fields
Info<< "Create neutro fields\n" << endl;

forAll(angularFlux_,energyI)
{
    angularFlux_[energyI].setSize(discreteDirections_);
    forAll(angularFlux_[energyI],directionI)
    {
        angularFlux_[energyI].set
        (
            directionI,
            new volScalarField
            (
                IOobject
                (
                    "angularFlux_"+Foam::name(energyI)+"_"+Foam::name(directionI),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                defaultFlux_
            )
        );
    }
}

forAll(flux_,fluxI)
{
    flux_.set
    (
        fluxI,
        new volScalarField
        (
            IOobject
            (
                "flux"+Foam::name(fluxI),
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            defaultFlux_
        )
    );
}


scalar pTot(power()); 
#include "normFluxesSN.H"
#include "updateScalarFluxes.H"

forAll(prec_,precI)
{
    prec_.set
    (
        precI,
        new volScalarField
        (
            IOobject
            (
                "prec"+Foam::name(precI),
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            defaultPrec_
        )
    );
}

forAll(prec_,precI)
{
    delayedNeutroSource_ += xs_.lambda()[precI]*prec_[precI];
}

forAll(facePhis_, directionI)
{
    facePhis_.set
    (
        directionI,
        new surfaceScalarField("",directionVersors_[directionI] & mesh_.Sf())
    );
}    




