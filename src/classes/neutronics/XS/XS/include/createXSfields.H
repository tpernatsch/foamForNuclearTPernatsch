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
                isReadWriteXS_ ? IOobject::READ_IF_PRESENT : IOobject::NO_READ,
                isReadWriteXS_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
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
                isReadWriteXS_ ? IOobject::READ_IF_PRESENT : IOobject::NO_READ,
                isReadWriteXS_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
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
                isReadWriteXS_ ? IOobject::READ_IF_PRESENT : IOobject::NO_READ,
                isReadWriteXS_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("", dimensionSet(1,1,-2,0,0,0,0), 1.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}

forAll(sigmaRemoval_, energyI)
{
    sigmaRemoval_.set
    (
        energyI,
        new volScalarField
        (
            IOobject
            (
                "sigmaRemoval"+Foam::name(energyI),
                mesh_.time().timeName(),
                mesh_,
                isReadWriteXS_ ? IOobject::READ_IF_PRESENT : IOobject::NO_READ,
                isReadWriteXS_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
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
                            isReadWriteXS_ ? IOobject::READ_IF_PRESENT : IOobject::NO_READ,
                            isReadWriteXS_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
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

// Create variable fields (e.g TFuel, rhoCool, ...) or read if already set
forAll(xsVariableNames_, varI)
{
    word variableName(xsVariableNames_[varI]);

    if (mesh_.foundObject<volScalarField>(variableName))
    {
        volScalarField& varField = const_cast<volScalarField&>
        (
            mesh_.lookupObject<volScalarField>(variableName)
        );

        xsVariableFields_.set(varI, &varField);
    }
    else
    {
        xsVariableFields_.set
        (
            varI,
            new volScalarField
            (
                IOobject
                (
                    variableName,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("", dimensionSet(0,0,0,0,0,0,0), 0.0),
                zeroGradientFvPatchScalarField::typeName
            )
        );
    }
}
