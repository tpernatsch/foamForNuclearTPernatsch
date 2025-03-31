IOobject alphaPhi10Header
(
    IOobject::groupName("alphaPhi0", alpha1.group()),
    dynMesh.time().timeName(),
    dynMesh,
    IOobject::READIFPRESENT,
    IOobject::AUTOWRITE
);

const bool alphaRestart =
    alphaPhi10Header.typeHeaderOk<surfaceScalarField>(true);

if (alphaRestart)
{
    Info << "Restarting alpha" << endl;
}

// MULES flux from previous time-step
surfaceScalarField alphaPhi10
(
    alphaPhi10Header,
    phi*fvc::interpolate(alpha1)
);

// MULES Correction
tmp<surfaceScalarField> talphaPhi1Corr0;

// MULES compressed flux is registered in case scalarTransport FO needs it.
surfaceScalarField alphaPhiUn
(
    IOobject
    (
        "alphaPhiUn",
        dynMesh.time().timeName(),
        dynMesh,
        IOobject::NOREAD,
        IOobject::NOWRITE
    ),
    dynMesh,
    dimensionedScalar(phi.dimensions(), Zero)
);
