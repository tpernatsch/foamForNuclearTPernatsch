volScalarField::Internal Sp
(
    IOobject
    (
        "Sp",
        runTime().timeName(),
        mesh_
    ),
    mesh_,
    dimensionedScalar(dgdt_.dimensions(), Zero)
);

volScalarField::Internal Su
(
    IOobject
    (
        "Su",
        runTime().timeName(),
        mesh_
    ),
    mesh_,
    dimensionedScalar(dgdt_.dimensions(), Zero)
);

forAll(dgdt_, celli)
{
    if (dgdt_[celli] > 0.0)
    {
        Sp[celli] -= dgdt_[celli]/max(1.0 - alpha1_[celli], 1e-4);
        Su[celli] += dgdt_[celli]/max(1.0 - alpha1_[celli], 1e-4);
    }
    else if (dgdt_[celli] < 0.0)
    {
        Sp[celli] += dgdt_[celli]/max(alpha1_[celli], 1e-4);
    }
}

volScalarField::Internal divU
(
    mesh_.moving()
  ? fvc::div(phiCN() + mesh_.phi())
  : fvc::div(phiCN())
);
