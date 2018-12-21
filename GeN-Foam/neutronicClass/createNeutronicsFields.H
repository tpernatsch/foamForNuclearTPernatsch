// Nuclear data fields
Info<< "Nuclear data fields\n" << endl;


forAll(IV, i)
{
    IV.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "IV"+Foam::name(i),
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

forAll(D, i)
{
    D.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "D"+Foam::name(i),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,//IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("", dimensionSet(0,1,0,0,0,0,0), 1.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}

forAll(nuSigmaEff, i)
{
    nuSigmaEff.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "nuSigmaEff"+Foam::name(i),
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

forAll(sigmaPow, i)
{
    sigmaPow.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "sigmaPow"+Foam::name(i),
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

forAll(sigmaDisapp, i)
{
    sigmaDisapp.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "sigmaDisapp"+Foam::name(i),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar("", dimensionSet(0,-1,0,0,0,0,0), 1.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}

forAll(sigmaFromTo,energyJ)
{
    sigmaFromTo[energyJ].resize(energyGroups); 

    forAll(sigmaFromTo[energyJ],energyI)
    {
        sigmaFromTo[energyJ].set
        (
            energyI,
            new volScalarField
            (
                IOobject
                (
                    "sigmaFromTo"+Foam::name(energyJ)+Foam::name(energyI),
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

forAll(chiPrompt, i)
{
    chiPrompt.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "chiPrompt"+Foam::name(i),
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

forAll(chiDelayed, i)
{
    chiDelayed.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "chiDelayed"+Foam::name(i),
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

forAll(Beta, i)
{
    Beta.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "Beta"+Foam::name(i),
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

forAll(lambda, i)
{
    lambda.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "lambda"+Foam::name(i),
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

forAll(discFactor,energyI)
{
    discFactor.set
    (
        energyI,
        new volScalarField
        (
            IOobject
            (
                "discFactor"+Foam::name(energyI),
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar("", dimensionSet(0,0,0,0,0,0,0), 1.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}

// Physical fields

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
    delayedNeutroSource += lambda[precI]*prec[precI];
}

#include "setNeutronicsConstants.H"
#include "setNeutronicsVariables.H"
#include "setPrecConst.H"





