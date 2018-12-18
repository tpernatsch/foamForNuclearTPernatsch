
volVectorField Disp
(
    IOobject
    (
        "Disp",
        runTime.timeName(),
        mechMesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mechMesh
);


volScalarField fuelDisp
(
    IOobject
    (
        "fuelDisp",
        runTime.timeName(),
        mechMesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mechMesh
);


volScalarField CRDisp
(
    IOobject
    (
        "CRDisp",
        runTime.timeName(),
        mechMesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mechMesh
);


volVectorField fuelDispVector
(
    IOobject
    (
        "fuelDispVector",
        runTime.timeName(),
        mechMesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    (fuelDisp*vector(0,0,1))
);

volVectorField meshDisp
(
    IOobject
    (
        "meshDisp",
        runTime.timeName(),
        mechMesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mechMesh,
    dimensionedVector("", dimensionSet(0,1,0,0,0,0,0), vector(0,0,0)),
    zeroGradientFvPatchScalarField::typeName
);

volVectorField meshDispOld
(
    IOobject
    (
        "meshDispOld",
        runTime.timeName(),
        mechMesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mechMesh,
    dimensionedVector("", dimensionSet(0,1,0,0,0,0,0), vector(0,0,0)),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField Tmech
(
    IOobject
    (
        "Tmech",
        runTime.timeName(),
        mechMesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mechMesh,
    TrefStructures,//important to impose correct BC (this value will remain in the soft region surrounding everything)
    zeroGradientFvPatchScalarField::typeName
);


Info<< "Calculating stress field sigmaD\n" << endl;
volSymmTensorField sigmaD
(
    IOobject
    (
        "sigmaD",
        runTime.timeName(),
        mechMesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    muMech*twoSymm(fvc::grad(Disp)) + lambdaMech*(I*tr(fvc::grad(Disp)))
);

Info<< "Calculating explicit part of div(sigma) divSigmaExp\n" << endl;
volVectorField divSigmaExp
(
    IOobject
    (
        "divSigmaExp",
        runTime.timeName(),
        mechMesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::div(sigmaD)
);

if (compactNormalStress)
{
    divSigmaExp -= fvc::laplacian(2*muMech + lambdaMech, Disp, "laplacian(DD,D)");
}
else
{
    divSigmaExp -= fvc::div((2*muMech + lambdaMech)*fvc::grad(Disp), "div(sigmaD)");
}

//create fields interpolated from the porous solver

volScalarField TavFuelMech
(
    IOobject
    (
        "TavFuelMech",
        runTime.timeName(),
        mechMesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mechMesh,
    dimensionedScalar("", dimensionSet(0,0,0,1,0,0,0), 1.0),
    zeroGradientFvPatchScalarField::typeName
);

