surfaceScalarField& phi1(fluid1_.phi());
surfaceScalarField& phi2(fluid2_.phi());

dimensionedScalar residualAlpha1(fluid1_.residualAlpha());
dimensionedScalar residualAlpha2(fluid2_.residualAlpha());

surfaceScalarField alphaf1(fvc::interpolate(alpha1));
surfaceScalarField alphaf2(fvc::interpolate(alpha2));
surfaceScalarField movingAlphaf(fvc::interpolate(movingAlpha_));

// Construct 1/A fields and limit to avoid problems for alpha = 0
volScalarField rAU1
(
    "rAU."+fluid1_.name(),
    1.0
   /(
        UEqn1.A()
      + max(residualAlpha1 - alpha1, 0.0)
       *rho1/mesh_.time().deltaT()
    )
);
volScalarField rAU2
(
    "rAU."+fluid2_.name(),
    1.0
   /(
        UEqn2.A()
      + max(residualAlpha2 - alpha2, 0.0)
       *rho2/mesh_.time().deltaT()
    )
);

// Interpolate 1/A to faces and limit
surfaceScalarField alpharAUf1
(
    fvc::interpolate(max(alpha1, residualAlpha1)*rAU1)
);
surfaceScalarField alpharAUf2
(
    fvc::interpolate(max(alpha2, residualAlpha2)*rAU2)
);

// --- Pressure corrector loop
while (pimple_.correct())
{
    this->correctContErrs();

    volScalarField rho(this->rho());

    // Correct p_rgh for consistency with p and the updated densities
    p_rgh_ = p_ - rho*gh_;

    // Correct fixed-flux BCs to be consistent with the velocity BCs
    MRF_.correctBoundaryFlux(U1, phi1);
    MRF_.correctBoundaryFlux(U2, phi2);

    // Construct H/A fields and limit to avoid problems for alpha = 0
    volVectorField HbyA1
    (
        "HbyA."+fluid1_.name(),
        U1
    );
    HbyA1 =
        rAU1
       *(
            UEqn1.H()
          + max(residualAlpha1 - alpha1, 0.0)
           *rho1*U1.oldTime()/mesh_.time().deltaT()
        );
    volVectorField HbyA2
    (
        "HbyA."+fluid2_.name(),
        U2
    );
    HbyA2 =
        rAU2
       *(
            UEqn2.H()
         +  max(residualAlpha2 - alpha2, 0.0)
           *rho2*U2.oldTime()/mesh_.time().deltaT()
        );

    // Construct fluxes due to buoyancy/gravity effects
    surfaceScalarField ghSnGradRho
    (
        "ghSnGradRho",
        ghf_*fvc::snGrad(rho)*mesh_.magSf()
    );
    surfaceScalarField phig1
    (
        alpharAUf1
       *(
           ghSnGradRho
         - alphaf2*fvc::interpolate(rho1 - rho2)*(g_ & mesh_.Sf())
        )
    );
    surfaceScalarField phig2
    (
        alpharAUf2
       *(
           ghSnGradRho
         - alphaf1*fvc::interpolate(rho2 - rho1)*(g_ & mesh_.Sf())
        )
    );

    // Predict phase fluxes
    surfaceScalarField phiHbyA1
    (
        "phiHbyA."+fluid1_.name(),
        fvc::flux(HbyA1)
      - phig1
    );
    surfaceScalarField phiHbyA2
    (
        "phiHbyA."+fluid2_.name(),
        fvc::flux(HbyA2)
      - phig2
    );

    volScalarField rAUKd1 = rAU1*Kd12;
    volScalarField rAUKd2 = rAU2*Kd12;

    // Face-drag coefficients
    surfaceScalarField rAUKdf1(fvc::interpolate(rAUKd1));
    surfaceScalarField rAUKdf2(fvc::interpolate(rAUKd2));

    // If NOT partial elimination
    if (!partialElimination_)
    {
        HbyA1 += rAUKd1*U2;
        HbyA2 += rAUKd2*U1;
        phiHbyA1 += rAUKdf1*MRF_.absolute(phi2);
        phiHbyA2 += rAUKdf2*MRF_.absolute(phi1);
    }
    
    // Construct the mean mixture predicted flux
    // including explicit drag contributions based on absolute fluxes
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        alphaf1*phiHbyA1
      + alphaf2*phiHbyA2
    );

    if (partialElimination_)
    {
        phiHbyA += 
            alphaf1*rAUKdf1*MRF_.absolute(phi2) +
            alphaf2*rAUKdf2*MRF_.absolute(phi1);
    }

    // If partial elimination
    //phiHbyA += alphafl*rAUKdl*MRF.absolute(phiv) + alphafv*rAUKdv*MRF.absolute(phil);

    // Interpolate mixture 1/A field
    surfaceScalarField rAUf
    (
        "rAUf",
        mag(alphaf1*alpharAUf1 + alphaf2*alpharAUf2)
    );

    // Update the fixedFluxPressure BCs to ensure flux consistency
    setSnGrad<fixedFluxPressureFvPatchScalarField>
    (
        p_rgh_.boundaryFieldRef(),
        (
            phiHbyA.boundaryField()
          - (
                alphaf1.boundaryField()*phi1.boundaryField()
              + alphaf2.boundaryField()*phi2.boundaryField()
            )
        )/(mesh_.magSf().boundaryField()*rAUf.boundaryField())
    );

    // Construct the compressibility parts of the pressure equation
    tmp<fvScalarMatrix> pEqnComp1;
    tmp<fvScalarMatrix> pEqnComp2;
    pEqnComp1 =
        (
            fluid1_.contErr() // +dmdt in contErr
          - fvc::Sp(fvc::ddt(alpha1) + fvc::div(fluid1_.alphaPhi()), rho1)
        )/rho1
      + (alpha1*fluid1_.thermo().psi()/rho1)*correction(fvm::ddt(p_rgh_));

    pEqnComp2 =
        (
            fluid2_.contErr() // -dmdt in contErr
          - fvc::Sp(fvc::ddt(alpha2) + fvc::div(fluid2_.alphaPhi()), rho2)
        )/rho2
      + (alpha2*fluid2_.thermo().psi()/rho2)*correction(fvm::ddt(p_rgh_));

    // Cache p prior to solve for density update
    volScalarField p_rgh_0(p_rgh_);

    // Iterate over the pressure equation to correct for non-orthogonality
    while (pimple_.correctNonOrthogonal())
    {
        fvScalarMatrix pEqnIncomp
        (
            fvc::div(phiHbyA)
          - fvm::laplacian(rAUf, p_rgh_)
        );

        fvScalarMatrix pEqn(pEqnComp1() + pEqnComp2() + pEqnIncomp);

        pEqn.setReference(pRefCell_, pRefValue_);

        solve
        (
            pEqn,
            mesh_.solver(p_rgh_.select(pimple_.finalInnerIter()))
        );

        if (pimple_.finalNonOrthogonalIter())
        {
            phi_ = phiHbyA + pEqnIncomp.flux();

            surfaceScalarField mSfGradp("mSfGradp", pEqnIncomp.flux()/rAUf);
            
            phi1 = phiHbyA1 + alpharAUf1*mSfGradp;
            phi2 = phiHbyA2 + alpharAUf2*mSfGradp;

            if (partialElimination_)
            {
                surfaceScalarField phir
                (
                    ((phi1 + rAUKdf1*phi2) - (phi2 + rAUKdf2*phi1))
                   /(1 - rAUKdf1*rAUKdf2)
                );

                phi1 = phi_/movingAlphaf + alphaf2*phir;
                phi2 = phi_/movingAlphaf - alphaf1*phir;
            }

            fluid1_.dgdt() = -pEqnComp1 & p_rgh_;
            fluid2_.dgdt() = -pEqnComp2 & p_rgh_;

            p_rgh_.relax();

            mSfGradp = pEqnIncomp.flux()/rAUf;

            volVectorField Us1 = 
                HbyA1 + fvc::reconstruct(alpharAUf1*mSfGradp - phig1);
            volVectorField Us2 = 
                HbyA2 + fvc::reconstruct(alpharAUf2*mSfGradp - phig2);

            if (partialElimination_)
            {
                //- It turns out that using:
                //  U = fvc::reconstruct(MRF_.absolute(phi_));
                //  was not correct, even though it solved single-phase
                //  oscillations. I will keep it in the single-phase solver.

                volVectorField U
                (
                    alpha1*(Us1 + rAUKd1*U2) + alpha2*(Us2 + rAUKd2*U1)
                );

                volVectorField Ur
                (
                    ((1 - rAUKd2)*Us1 - (1 - rAUKd1)*Us2)/(1 - rAUKd1*rAUKd2)
                );

                U1 = U/movingAlpha_ + alpha2*Ur;
                U2 = U/movingAlpha_ - alpha1*Ur;
            }
            else
            {
                U1 = Us1;
                U2 = Us2;
            }

            U1.correctBoundaryConditions();
            U2.correctBoundaryConditions();
            fvOptions_.correct(U1);
            fvOptions_.correct(U2);
            fluid1_.magU() = mag(U1);
            fluid2_.magU() = mag(U2);
        }
    }
    
    // Update and limit the static pressure
    p_ = max(p_rgh_ + rho*gh_, pMin_);

    // Limit p_rgh
    p_rgh_ = p_ - rho*gh_;

    // Update densities from change in p_rgh
    rho1 += fluid1_.thermo().psi()*(p_rgh_ - p_rgh_0);
    rho2 += fluid2_.thermo().psi()*(p_rgh_ - p_rgh_0);

    // Correct p_rgh for consistency with p and the updated densities
    rho = rho1*alpha1 + rho2*alpha2;
    p_rgh_ = p_ - rho*gh_;
    p_rgh_.correctBoundaryConditions();
}

{
    this->correctContErrs();
    volScalarField contErrRel1(fluid1_.contErr()/rho1);
    volScalarField contErrRel2(fluid2_.contErr()/rho2);

    Info    << "Continuity error " << fluid1_.name() << " (avg min max) ="
            << " " << contErrRel1.weightedAverage(mesh_.V()).value()
            << " " << min(contErrRel1).value()
            << " " << max(contErrRel1).value()
            << " 1/s"  << endl;
    Info    << "Continuity error " << fluid2_.name() << " (avg min max) ="
            << " " << contErrRel2.weightedAverage(mesh_.V()).value()
            << " " << min(contErrRel2).value()
            << " " << max(contErrRel2).value() 
            << " 1/s" << endl;
}
