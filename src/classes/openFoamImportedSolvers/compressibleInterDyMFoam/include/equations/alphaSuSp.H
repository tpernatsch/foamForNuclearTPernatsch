volScalarField::Internal Sp
(
    IOobject
    (
        "Sp",
        runTime.timeName(),
        dynMesh
    ),
    dynMesh,
    dimensionedScalar(dgdt.dimensions(), Zero)
);

volScalarField::Internal Su
(
    IOobject
    (
        "Su",
        runTime.timeName(),
        dynMesh
    ),
    dynMesh,
    dimensionedScalar(dgdt.dimensions(), Zero)
);

forAll(dgdt, celli)
{
    if (dgdt[celli] > 0.0 && alpha1[celli] > 0.0)
    {
        Sp[celli] -= dgdt[celli]*alpha1[celli];
        Su[celli] += dgdt[celli]*alpha1[celli];
    }
    else if (dgdt[celli] < 0.0 && alpha1[celli] < 1.0)
    {
        Sp[celli] += dgdt[celli]*(1.0 - alpha1[celli]);
    }
}

volScalarField::Internal divU
(
    dynMesh.moving()
  ? fvc::div(phiCN() + dynMesh.phi())
  : fvc::div(phiCN())
);
