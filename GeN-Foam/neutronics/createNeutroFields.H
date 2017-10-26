// Nuclear data fields
Info<< "Nuclear data fields\n" << endl;

PtrList<volScalarField> IV(energyGroups);
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
                runTime.timeName(),
                neutroMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            neutroMesh,
            dimensionedScalar("", dimensionSet(0,-1,1,0,0,0,0), 1.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}

PtrList<volScalarField> D(energyGroups);
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
                runTime.timeName(),
                neutroMesh,
                IOobject::NO_READ,//IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            neutroMesh,
            dimensionedScalar("", dimensionSet(0,1,0,0,0,0,0), 1.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}

volScalarField Dalbedo
(

    IOobject
    (
        "Dalbedo",
        runTime.timeName(),
        neutroMesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    neutroMesh,
    dimensionedScalar("", dimensionSet(0,1,0,0,0,0,0), 1.0),
    zeroGradientFvPatchScalarField::typeName

);

PtrList<volScalarField> nuSigmaEff(energyGroups);
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
                runTime.timeName(),
                neutroMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            neutroMesh,
            dimensionedScalar("", dimensionSet(0,-1,0,0,0,0,0), 1.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}

PtrList<volScalarField> sigmaPow(energyGroups);
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
                runTime.timeName(),
                neutroMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            neutroMesh,
            dimensionedScalar("", dimensionSet(1,1,-2,0,0,0,0), 1.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}

PtrList<volScalarField> sigmaDisapp(energyGroups);
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
                runTime.timeName(),
                neutroMesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            neutroMesh,
            dimensionedScalar("", dimensionSet(0,-1,0,0,0,0,0), 1.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}

List<PtrList<volScalarField > > sigmaFromTo(energyGroups);

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
                    runTime.timeName(),
                    neutroMesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                neutroMesh,
                dimensionedScalar("", dimensionSet(0,-1,0,0,0,0,0), 1.0),
                zeroGradientFvPatchScalarField::typeName
            )
        );
    }

}

PtrList<volScalarField> chiPrompt(energyGroups);
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
                runTime.timeName(),
                neutroMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            neutroMesh,
            dimensionedScalar("", dimensionSet(0,0,0,0,0,0,0), 1.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}

PtrList<volScalarField> chiDelayed(energyGroups);
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
                runTime.timeName(),
                neutroMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            neutroMesh,
            dimensionedScalar("", dimensionSet(0,0,0,0,0,0,0), 1.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}

PtrList<volScalarField> Beta(precGroups);
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
                runTime.timeName(),
                neutroMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            neutroMesh,
            dimensionedScalar("", dimensionSet(0,0,0,0,0,0,0), 1.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}

volScalarField BetaTot
(
    IOobject
    (
        "BetaTot",
        runTime.timeName(),
        neutroMesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    neutroMesh,
    dimensionedScalar("", dimensionSet(0,0,0,0,0,0,0), 1.0),
    zeroGradientFvPatchScalarField::typeName
);

PtrList<volScalarField> lambda(precGroups);
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
                runTime.timeName(),
                neutroMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            neutroMesh,
            dimensionedScalar("", dimensionSet(0,0,-1,0,0,0,0), 1.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}

volScalarField fuelFraction
(
    IOobject
    (
        "fuelFraction",
        runTime.timeName(),
        neutroMesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    neutroMesh,
    dimensionedScalar("", dimensionSet(0,0,0,0,0,0,0), 1.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField ScNo
(
    IOobject
    (
        "ScNo",
        runTime.timeName(),
        neutroMesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    neutroMesh,
    dimensionedScalar("", dimensionSet(0,0,0,0,0,0,0), 1.0),
    zeroGradientFvPatchScalarField::typeName
);

PtrList<volScalarField > discFactor(energyGroups);
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
                runTime.timeName(),
                neutroMesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            neutroMesh,
            dimensionedScalar("", dimensionSet(0,0,0,0,0,0,0), 1.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}


// Physical fields

PtrList<volScalarField> flux(energyGroups);
PtrList<volScalarField> fluxStar(energyGroups);//flux*discFactor
PtrList<volScalarField> fluxStar2(energyGroups);//flux*discFactor
PtrList<volScalarField> prec(precGroups);


volScalarField fluxStarAlbedo
(
    IOobject
    (
        "fluxStarAlbedo",
        runTime.timeName(),
        neutroMesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
            neutroMesh,
    dimensionedScalar("", dimensionSet(0,-2,-1,0,0,0,0), 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField defaultFlux
(
    IOobject
    (
        "defaultFlux",
        runTime.timeName(),
        neutroMesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    neutroMesh
);

volScalarField defaultFlux2
(
    IOobject
    (
        "defaultFlux2",
        runTime.timeName(),
        neutroMesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    neutroMesh
);

volScalarField defaultPrec
(
    IOobject
    (
        "defaultPrec",
        runTime.timeName(),
        neutroMesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
            neutroMesh,
    dimensionedScalar("", dimensionSet(0,-3,0,0,0,0,0), 0.0),
    zeroGradientFvPatchScalarField::typeName
);

forAll(flux,fluxI)
{
    Info<< "    Adding to flux\n" << endl;
    flux.set
    (
        fluxI,
        new volScalarField
        (
            IOobject
            (
                "flux"+Foam::name(fluxI),
                runTime.timeName(),
                neutroMesh,
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
                    runTime.timeName(),
                    neutroMesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                defaultFlux
            )
        );
}

List<PtrList<volScalarField > > base_fluxStar(energyGroups);

forAll(fluxStar,energyI)
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

forAll(fluxStar2,fluxI)
{
        fluxStar2.set
        (
            fluxI,
            new volScalarField
            (
                    IOobject
                    (
                        "fluxStar2"+Foam::name(fluxI),
                        runTime.timeName(),
                        neutroMesh,
                        IOobject::READ_IF_PRESENT,
                        IOobject::AUTO_WRITE
                    ),
                    defaultFlux2
                )
        );
}

volScalarField oneGroupFlux
(
    IOobject
    (
        "oneGroupFlux",
        runTime.timeName(),
        neutroMesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    defaultFlux
);
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
                runTime.timeName(),
                neutroMesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            defaultPrec
        )
    );
}

volScalarField neutroSource
(
    IOobject
    (
        "neutroSource",
        runTime.timeName(),
        neutroMesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    neutroMesh,
    dimensionedScalar("", dimensionSet(0,-3,-1,0,0,0,0), 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField delayedNeutroSource
(
    IOobject
    (
        "delayedNeutroSource",
        runTime.timeName(),
        neutroMesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    neutroMesh,
    dimensionedScalar("", dimensionSet(0,-3,-1,0,0,0,0), 0.0),
    zeroGradientFvPatchScalarField::typeName
);
forAll(prec,precI)
{
    delayedNeutroSource += lambda[precI]*prec[precI];
}

volScalarField scatteringSourceExtra
(
    IOobject
    (
        "scatteringSourceExtra",
        runTime.timeName(),
        neutroMesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    neutroMesh,
    dimensionedScalar("", dimensionSet(0,-3,-1,0,0,0,0), 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volVectorField UNeutro
(
    IOobject
    (
        "UNeutro",
        runTime.timeName(),
        neutroMesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    neutroMesh,
    dimensionedVector("", dimensionSet(0,1,-1,0,0,0,0), vector(0.0,0.0,0.0)),
    zeroGradientFvPatchVectorField::typeName
);

surfaceScalarField phiNeutro
(
    IOobject
    (
        "phiNeutro",
        runTime.timeName(),
        neutroMesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
            linearInterpolate(UNeutro) & neutroMesh.Sf()
);

volScalarField porosityNeutro
(
    IOobject
    (
        "porosityNeutro",
        runTime.timeName(),
        neutroMesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    neutroMesh,
    dimensionedScalar("", dimensionSet(0,0,0,0,0,0,0),1.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField volFuelPowerNeutro
(
    IOobject
    (
        "volFuelPowerNeutro",
        runTime.timeName(),
        neutroMesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    neutroMesh,
    dimensionedScalar("", dimensionSet(1,-1,-3,0,0,0,0), 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField TfuelNeutro
(
    IOobject
    (
        "TfuelNeutro",
        runTime.timeName(),
        neutroMesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    neutroMesh,
    dimensionedScalar("", dimensionSet(0,0,0,1,0,0,0),TfuelRef),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField TcladNeutro
(
    IOobject
    (
        "TcladNeutro",
        runTime.timeName(),
        neutroMesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    neutroMesh,
    dimensionedScalar("", dimensionSet(0,0,0,1,0,0,0),TcladRef),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField radExp
(
    IOobject
    (
        "radExp",
        runTime.timeName(),
        neutroMesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    neutroMesh,
    dimensionedScalar("", dimensionSet(0,0,0,0,0,0,0), 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField axExp
(
    IOobject
    (
        "axExp",
        runTime.timeName(),
        neutroMesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    neutroMesh,
    dimensionedScalar("", dimensionSet(0,0,0,0,0,0,0), 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField rhoCoolNeutro
(
    IOobject
    (
        "rhoCoolNeutro",
        runTime.timeName(),
        neutroMesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    neutroMesh,
    dimensionedScalar("", dimensionSet(1,-3,0,0,0,0,0), rhoCoolRef),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField TCoolNeutro
(
    IOobject
    (
        "TCoolNeutro",
        runTime.timeName(),
        neutroMesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    neutroMesh,
    dimensionedScalar("", dimensionSet(0,0,0,1,0,0,0),TCoolRef),
    zeroGradientFvPatchScalarField::typeName
);

//displacement field for the mesh (should come from mechanical solver)
volVectorField DNeutro
(
    IOobject
    (
        "DNeutro",
        runTime.timeName(),
        neutroMesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    neutroMesh,
    dimensionedVector("d_zero", dimensionSet(0,1,0,0,0,0,0), vector(0,0,0)),
    zeroGradientFvPatchScalarField::typeName
);

#include "setNeutronicsConstants.H"
#include "setNeutronicsVariables.H"
#include "setPrecConst.H"

//for aitken
List<scalar> integralFlux0(energyGroups);
List<scalar> integralFlux1(energyGroups);
List<scalar> integralFlux2(energyGroups);

//for POD acc
scalar neutroResidual0 =0.0;
scalar neutroResidual1 =0.0;
scalar neutroResidual2 =0.0;



