// Nuclear data fields
Info<< "Creating nuclear data fields\n" << endl;

forAll(IV_, energyI)
{
    IV_.set
    (
        energyI,
        new volScalarField
        (
            IOobject
            (
                "IV"+Foam::name(energyI),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("", dimensionSet(0,-1,1,0,0,0,0), 1.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}

forAll(D_, energyI)
{
    D_.set
    (
        energyI,
        new volScalarField
        (
            IOobject
            (
                "D"+Foam::name(energyI),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("", dimensionSet(0,1,0,0,0,0,0), 1.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}

forAll(nuSigmaEff_, energyI)
{
    nuSigmaEff_.set
    (
        energyI,
        new volScalarField
        (
            IOobject
            (
                "nuSigmaEff"+Foam::name(energyI),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("", dimensionSet(0,-1,0,0,0,0,0), 1.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}

forAll(sigmaPow_, energyI)
{
    sigmaPow_.set
    (
        energyI,
        new volScalarField
        (
            IOobject
            (
                "sigmaPow"+Foam::name(energyI),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("", dimensionSet(1,1,-2,0,0,0,0), 1.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}

forAll(sigmaDisapp_, energyI)
{
    sigmaDisapp_.set
    (
        energyI,
        new volScalarField
        (
            IOobject
            (
                "sigmaDisapp"+Foam::name(energyI),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("", dimensionSet(0,-1,0,0,0,0,0), 1.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}

if (!isLowMemory_)
{
    forAll(sigmaFromTo_, momentI)
    {
        sigmaFromTo_[momentI].resize(energyGroups_);
        forAll(sigmaFromTo_[momentI], energyJ)
        {
            sigmaFromTo_[momentI][energyJ].resize(energyGroups_); 
            forAll(sigmaFromTo_[momentI][energyJ], energyI)
            {
                sigmaFromTo_[momentI][energyJ].set
                (
                    energyI,
                    new volScalarField
                    (
                        IOobject
                        (
                            "sigmaFromToP"+Foam::name(momentI)+"_"+Foam::name(energyJ)+Foam::name(energyI),
                            mesh_.time().timeName(),
                            mesh_,
                            IOobject::NO_READ,
                            IOobject::NO_WRITE
                        ),
                        mesh_,
                        dimensionedScalar("", dimensionSet(0,-1,0,0,0,0,0), 1.0),
                        zeroGradientFvPatchScalarField::typeName
                    )
                );
            }
        } 
    }
}

forAll(chiPrompt_, energyI)
{
    chiPrompt_.set
    (
        energyI,
        new volScalarField
        (
            IOobject
            (
                "chiPrompt"+Foam::name(energyI),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("", dimensionSet(0,0,0,0,0,0,0), 1.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}

forAll(chiDelayed_, energyI)
{
    chiDelayed_.set
    (
        energyI,
        new volScalarField
        (
            IOobject
            (
                "chiDelayed"+Foam::name(energyI),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("", dimensionSet(0,0,0,0,0,0,0), 1.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}

forAll(Beta_, precI)
{
    Beta_.set
    (
        precI,
        new volScalarField
        (
            IOobject
            (
                "Beta"+Foam::name(precI),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("", dimensionSet(0,0,0,0,0,0,0), 1.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}

forAll(lambda_, precI)
{
    lambda_.set
    (
        precI,
        new volScalarField
        (
            IOobject
            (
                "lambda"+Foam::name(precI),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("", dimensionSet(0,0,-1,0,0,0,0), 1.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}

forAll(discFactor_, energyI)
{
    discFactor_.set
    (
        energyI,
        new volScalarField
        (
            IOobject
            (
                "discFactor"+Foam::name(energyI),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("", dimensionSet(0,0,0,0,0,0,0), 1.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}
