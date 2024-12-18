IOobject alphaPhi10Header
(
    IOobject::groupName("alphaPhi0", alpha1_.group()),
    mesh_.time().timeName(),
    mesh_,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
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
    phi_*fvc::interpolate(alpha1_)
);

// MULES Correction
tmp<surfaceScalarField> talphaPhi1Corr0;

// MULES compressed flux is registered in case scalarTransport FO needs it.
surfaceScalarField alphaPhiUn
(
    IOobject
    (
        "alphaPhiUn",
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh_,
    dimensionedScalar(phi_.dimensions(), Zero)
);
