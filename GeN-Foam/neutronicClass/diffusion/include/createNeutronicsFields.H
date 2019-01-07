// Nuclear data fields
Info<< "Nuclear data fields\n" << endl;

forAll(flux,fluxI)
{
    flux.set
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
            defaultFlux
        )
    );
}

forAll(fluxStar,fluxI)
{
        fluxStar.set
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
                defaultFlux
            )
        );
}

oneGroupFlux *= 0.0;
forAll(flux,energyI)
{
    oneGroupFlux+=flux[energyI];
}

forAll(prec,precI)
{
    prec.set
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
            defaultPrec
        )
    );
}

forAll(prec,precI)
{
    delayedNeutroSource += xs.lambda()[precI]*prec[precI];
}






