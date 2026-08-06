{
    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U_, p_rgh_));

    surfaceScalarField phig(-rAUf*ghf_*fvc::snGrad(rhok_)*mesh_.magSf());

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + MRF_.zeroFilter(rAUf*fvc::ddtCorr(U_, phi_))
      + phig
    );

    MRF_.makeRelative(phiHbyA);

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh_, U_, phiHbyA, rAUf, MRF_);

    while (pimple_.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUf, p_rgh_) == fvc::div(phiHbyA)
        );

        p_rghEqn.setReference(pRefCell_, getRefCellValue(p_rgh_, pRefCell_));

        p_rghEqn.solve(p_rgh_.select(pimple_.finalInnerIter()));

        if (pimple_.finalNonOrthogonalIter())
        {
            // Calculate the conservative fluxes
            phi_ = phiHbyA - p_rghEqn.flux();

            // Explicitly relax pressure for momentum corrector
            p_rgh_.relax();

            // Correct the momentum source with the pressure gradient flux
            // calculated from the relaxed pressure
            U_ = HbyA + rAU*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
            U_.correctBoundaryConditions();
            fvOptions_.correct(U_);

            // Correct Uf if the mesh is moving
            fvc::correctUf(Uf_, U_, phi_);

            // Make the fluxes relative to the mesh motion
            fvc::makeRelative(phi_, U_);
        }
    }

    {
        volScalarField contErr(fvc::div(phi_));

        scalar sumLocalContErr = mesh_.time().deltaTValue()*
            mag(contErr)().weightedAverage(mesh_.V()).value();

        scalar globalContErr = mesh_.time().deltaTValue()*
            contErr.weightedAverage(mesh_.V()).value();
        cumulativeContErr_ += globalContErr;

        Info<< "time step continuity errors : sum local = " << sumLocalContErr
            << ", global = " << globalContErr
            << ", cumulative = " << cumulativeContErr_
            << endl;
    }

    p_ = p_rgh_ + rhok_*gh_;

    if (p_rgh_.needReference())
    {
        p_ += dimensionedScalar
        (
            "p",
            p_.dimensions(),
            pRefValue_ - getRefCellValue(p_, pRefCell_)
        );
        p_rgh_ = p_ - rhok_*gh_;
    }
}
