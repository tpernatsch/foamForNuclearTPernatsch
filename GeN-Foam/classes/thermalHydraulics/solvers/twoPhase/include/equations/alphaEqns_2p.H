{
    const dictionary& alphaDict(mesh_.solverDict("alpha"));
    volScalarField& alpha1(fluid1_);
    volScalarField& alpha2(fluid2_);
    volScalarField& alphas(structure_);

    //- This is a lambda function (standard C++ stuff) that adjusts the
    //  boundary field of a given volumetric flux alphaPhi so that it is
    //  consistent with the BCs for the fields (alpha and U, thus phi) from
    //  which it is derived. It is basically a physically-motivated
    //  correctBoundaryConditions() for alphaPhi
    auto correctInflowOutflow = []
    (
        const volScalarField& alpha,
        const surfaceScalarField& phi,
        surfaceScalarField& alphaPhi
    )
    {
        const volScalarField::Boundary& alphaBf = alpha.boundaryField();
        const surfaceScalarField::Boundary& phiBf = phi.boundaryField();
        surfaceScalarField::Boundary& alphaPhiBf = alphaPhi.boundaryFieldRef();

        forAll(alphaPhiBf, patchi)
        {
            fvsPatchScalarField& alphaPhip = alphaPhiBf[patchi];

            if (!alphaPhip.coupled())
            {
                alphaPhip = phiBf[patchi]*alphaBf[patchi];
            }
        }
    };

    if (alphaEqnsSolver_ == alphaEqnsSolver::implicitUpwind)
    {
        #include "solveAlphaImplicitUpwind.H"
    }
    else 
    {
        #include "solveAlphaMULES.H"
    }

    //- Update normalized alpha fields (only useful for regime map 
    //  parametrization, if selected as regime map parameter)
    fluid1_.normalized() = fluid1_/movingAlpha_;
    fluid2_.normalized() = fluid2_/movingAlpha_;

    //- Correct flow qualities, Lockhart-Martinelli parameter
    fluid1_.correctTwoPhaseFields(fluid2_);
    fluid2_.correctTwoPhaseFields(fluid1_);

    //- Correct continuity errors
    correctContErrs();

    //- Infos
    Info<< alpha1.name() << " (avg min max) ="
        << " " << alpha1.weightedAverage(mesh_.V()).value()
        << " " << min(alpha1).value()
        << " " << max(alpha1).value()
        << endl;

    Info<< alpha2.name() << " (avg min max) ="
        << " " << alpha2.weightedAverage(mesh_.V()).value()
        << " " << min(alpha2).value()
        << " " << max(alpha2).value()
        << endl;

    volScalarField alphaSum(alpha1 + alpha2 + alphas);
    Info<< "sum(alpha)" << " (avg min max) ="
        << " " << alphaSum.weightedAverage(mesh_.V()).value()
        << " " << min(alphaSum).value()
        << " " << max(alphaSum).value()
        << endl;
}
