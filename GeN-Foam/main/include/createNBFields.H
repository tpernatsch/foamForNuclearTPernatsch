volScalarField TNB
(
    IOobject
    (
        "TNB",
        runTime.timeName(),
        fluidMeshNB,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fluidMeshNB,
    dimensionedScalar("", dimTemperature, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField TavFuelNB
(
    IOobject
    (
        "TavFuelNB",
        runTime.timeName(),
        fluidMeshNB,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fluidMeshNB,
    dimensionedScalar("",  dimTemperature, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField TavCladNB
(
    IOobject
    (
        "TavCladNB",
        runTime.timeName(),
        fluidMeshNB,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fluidMeshNB,
    dimensionedScalar("",  dimTemperature, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField TstructuresNB
(
    IOobject
    (
        "TstructuresNB",
        runTime.timeName(),
        fluidMeshNB,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fluidMeshNB,
    dimensionedScalar("", dimTemperature, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField porosityNB
(
    IOobject
    (
        "porosityNB",
        runTime.timeName(),
        fluidMeshNB,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fluidMeshNB,
    dimensionedScalar("", dimless, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField rhoNB
(
    IOobject
    (
        "rhoNB",
        runTime.timeName(),
        fluidMeshNB,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fluidMeshNB,
    dimensionedScalar("", dimMass/dimVol, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volVectorField UNB
(
    IOobject
    (
        "UNB",
        runTime.timeName(),
        fluidMesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fluidMeshNB,
    dimensionedVector("", dimVelocity, vector::zero),
    zeroGradientFvPatchScalarField::typeName
);
