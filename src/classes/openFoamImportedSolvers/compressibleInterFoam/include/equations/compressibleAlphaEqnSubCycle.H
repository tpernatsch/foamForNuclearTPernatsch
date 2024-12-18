tmp<surfaceScalarField> talphaPhi1(alphaPhi10);

if (nAlphaSubCycles > 1)
{
    dimensionedScalar totalDeltaT = runTime.deltaT();

    talphaPhi1.reset
    (
        new surfaceScalarField
        (
            mesh_.newIOobject("alphaPhi1"),
            mesh_,
            dimensionedScalar(alphaPhi10.dimensions(), Zero)
        )
    );

    surfaceScalarField rhoPhiSum
    (
        mesh_.newIOobject("rhoPhiSum"),
        mesh_,
        dimensionedScalar(rhoPhi_.dimensions(), Zero)
    );

    tmp<volScalarField> trSubDeltaT;

    if (LTS)
    {
        trSubDeltaT =
            fv::localEulerDdt::localRSubDeltaT(mesh_, nAlphaSubCycles);
    }

    for
    (
        subCycle<volScalarField> alphaSubCycle(alpha1_, nAlphaSubCycles);
        !(++alphaSubCycle).end();
    )
    {
        #include "alphaEqn.H"
        talphaPhi1.ref() += (runTime.deltaT()/totalDeltaT)*alphaPhi10;
        rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi_;
    }

    rhoPhi_ = rhoPhiSum;
}
else
{
    #include "alphaEqn.H"
}

rho_ == alpha1_*rho1_ + alpha2_*rho2_;

const surfaceScalarField& alphaPhi1 = talphaPhi1();
surfaceScalarField alphaPhi2("alphaPhi2", phi_ - alphaPhi1);

volScalarField::Internal contErr
(
    (
        fvc::ddt(rho_) + fvc::div(rhoPhi_)
      - (fvOptions_(alpha1_, mixture_.thermo1().rho())&rho1_)
      - (fvOptions_(alpha2_, mixture_.thermo2().rho())&rho2_)
    )()
);
