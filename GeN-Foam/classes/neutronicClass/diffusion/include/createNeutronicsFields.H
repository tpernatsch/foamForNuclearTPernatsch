// Nuclear data fields
Info<< "Create neutro fields\n" << endl;

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

forAll(fluxStar_,fluxI)
{
        fluxStar_.set
        (
            fluxI,
            new volScalarField
            (
                IOobject
                (
                    "fluxStar"+Foam::name(fluxI),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                defaultFlux_
            )
        );
}

oneGroupFlux_ *= 0.0;
forAll(flux_,energyI)
{
    oneGroupFlux_+=flux_[energyI];
}

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






