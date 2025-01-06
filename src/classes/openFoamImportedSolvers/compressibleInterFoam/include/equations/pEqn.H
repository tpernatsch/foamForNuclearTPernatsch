{
    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U_, p_rgh_));
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + MRF_.zeroFilter(fvc::interpolate(rho_*rAU)*fvc::ddtCorr(U_, phi_))
    );
    MRF_.makeRelative(phiHbyA);

    surfaceScalarField phig
    (
        (
            mixture_.surfaceTensionForce()
          - ghf_*fvc::snGrad(rho_)
        )*rAUf*mesh_.magSf()
    );

    phiHbyA += phig;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh_, U_, phiHbyA, rAUf, MRF_);

    tmp<fvScalarMatrix> p_rghEqnComp1;
    tmp<fvScalarMatrix> p_rghEqnComp2;

    if (pimple_.transonic())
    {
        #include "rhofs.H"

        surfaceScalarField phid1("phid1", fvc::interpolate(psi1_)*phi_);
        surfaceScalarField phid2("phid2", fvc::interpolate(psi2_)*phi_);

        p_rghEqnComp1 =
            pos(alpha1_)
           *(
                (
                    fvc::ddt(alpha1_, rho1_) + fvc::div(alphaPhi1*rho1f)
                  - (fvOptions_(alpha1_, mixture_.thermo1().rho())&rho1_)
                )/rho1_
              - fvc::ddt(alpha1_) - fvc::div(alphaPhi1)
              + (alpha1_/rho1_)
               *correction
                (
                    psi1_*fvm::ddt(p_rgh_)
                  + fvm::div(phid1, p_rgh_) - fvm::Sp(fvc::div(phid1), p_rgh_)
                )
            );
        p_rghEqnComp1.ref().relax();

        p_rghEqnComp2 =
            pos(alpha2_)
           *(
               (
                   fvc::ddt(alpha2_, rho2_) + fvc::div(alphaPhi2*rho2f)
                 - (fvOptions_(alpha2_, mixture_.thermo2().rho())&rho2_)
               )/rho2_
             - fvc::ddt(alpha2_) - fvc::div(alphaPhi2)
             + (alpha2_/rho2_)
              *correction
               (
                   psi2_*fvm::ddt(p_rgh_)
                 + fvm::div(phid2, p_rgh_) - fvm::Sp(fvc::div(phid2), p_rgh_)
               )
           );
        p_rghEqnComp2.ref().relax();
    }
    else
    {
        #include "rhofs.H"

        p_rghEqnComp1 =
            pos(alpha1_)
           *(
                (
                    fvc::ddt(alpha1_, rho1_) + fvc::div(alphaPhi1*rho1f)
                  - (fvOptions_(alpha1_, mixture_.thermo1().rho())&rho1_)
                )/rho1_
              - fvc::ddt(alpha1_) - fvc::div(alphaPhi1)
              + (alpha1_*psi1_/rho1_)*correction(fvm::ddt(p_rgh_))
            );

        p_rghEqnComp2 =
            pos(alpha2_)
           *(
               (
                   fvc::ddt(alpha2_, rho2_) + fvc::div(alphaPhi2*rho2f)
                 - (fvOptions_(alpha2_, mixture_.thermo2().rho())&rho2_)
               )/rho2_
             - fvc::ddt(alpha2_) - fvc::div(alphaPhi2)
             + (alpha2_*psi2_/rho2_)*correction(fvm::ddt(p_rgh_))
            );
    }

    // Cache p_rgh_ prior to solve for density update
    volScalarField p_rgh_0(p_rgh_);

    while (pimple_.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqnIncomp
        (
            fvc::div(phiHbyA)
          - fvm::laplacian(rAUf, p_rgh_)
        );

        solve
        (
            p_rghEqnComp1() + p_rghEqnComp2() + p_rghEqnIncomp,
            p_rgh_.select(pimple_.finalInnerIter())
        );

        if (pimple_.finalNonOrthogonalIter())
        {
            p_ = max(p_rgh_ + (alpha1_*rho1_ + alpha2_*rho2_)*gh_, pMin_);
            p_rgh_ = p_ - (alpha1_*rho1_ + alpha2_*rho2_)*gh_;

            dgdt_ =
            (
                alpha1_*(p_rghEqnComp2 & p_rgh_)
              - alpha2_*(p_rghEqnComp1 & p_rgh_)
            );

            phi_ = phiHbyA + p_rghEqnIncomp.flux();

            U_ = HbyA
              + rAU*fvc::reconstruct((phig + p_rghEqnIncomp.flux())/rAUf);
            U_.correctBoundaryConditions();
            fvOptions_.correct(U_);
        }
    }

    // Update densities from change in p_rgh_
    mixture_.thermo1().correctRho(psi1_*(p_rgh_ - p_rgh_0));
    mixture_.thermo2().correctRho(psi2_*(p_rgh_ - p_rgh_0));

    rho_ = alpha1_*rho1_ + alpha2_*rho2_;

    // Correct p_rgh_ for consistency with p_ and the updated densities
    p_rgh_ = p_ - rho_*gh_;
    p_rgh_.correctBoundaryConditions();

    K_ = 0.5*magSqr(U_);
}
