surfaceScalarField& phi1(fluid1_.phi());
surfaceScalarField& phi2(fluid2_.phi());

dimensionedScalar residualAlpha1(fluid1_.residualAlpha());
dimensionedScalar residualAlpha2(fluid2_.residualAlpha());

surfaceScalarField alphaf1(fvc::interpolate(alpha1));
surfaceScalarField alphaf2(fvc::interpolate(alpha2));
surfaceScalarField movingAlphaf(fvc::interpolate(movingAlpha_));

surfaceScalarField alphaRhof10
(
    "alphaRhof10",
    fvc::interpolate
    (
        max(alpha1.oldTime(), residualAlpha1)
       *rho1.oldTime()
    )
);
surfaceScalarField alphaRhof20
(
    "alphaRhof20",
    fvc::interpolate
    (
        max(alpha2.oldTime(), residualAlpha2)
       *rho2.oldTime()
    )
);

// Construct 1/A fields and limit to avoid problems for alpha = 0
surfaceScalarField rAUf1
(
    "rAUf."+fluid1_.name(),
    1.0
   /(
            fvc::interpolate(UEqn1.A())
       +    alphaRhof10/mesh_.time().deltaT()
    )
);
surfaceScalarField rAUf2
(
    "rAUf."+fluid2_.name(),
    1.0
   /(
            fvc::interpolate(UEqn2.A())
       +    alphaRhof20/mesh_.time().deltaT()
    )
);

// Interpolate 1/A to faces and limit
surfaceScalarField alpharAUf1
(
    max(alphaf1, residualAlpha1)*rAUf1
);
surfaceScalarField alpharAUf2
(
    max(alphaf2, residualAlpha2)*rAUf2
);

surfaceScalarField Kdf12(fvc::interpolate(Kd12));

while (pimple_.correct())
{
    this->correctContErrs();

    volScalarField rho(this->rho());

    // Correct p_rgh for consistency with p and the updated densities
    p_rgh_ = p_ - rho*gh_;

    surfaceScalarField rhof1(fvc::interpolate(rho1));
    surfaceScalarField rhof2(fvc::interpolate(rho2));

    // Correct fixed-flux BCs to be consistent with the velocity BCs
    MRF_.correctBoundaryFlux(U1, phi1);
    MRF_.correctBoundaryFlux(U2, phi2);

    //- Compute buoyancy fluxes
    surfaceScalarField ghSnGradRho
    (
        "ghSnGradRho",
        ghf_*fvc::snGrad(rho)*mesh_.magSf()
    );
    surfaceScalarField phig1
    (
        "phig."+fluid1_.name(),
        alpharAUf1
       *(
            ghSnGradRho
          - alphaf2*(rhof1 - rhof2)*(g_ & mesh_.Sf())
        )
    );
    surfaceScalarField phig2
    (
        "phig."+fluid2_.name(),
        alpharAUf2
       *(
            ghSnGradRho
          - alphaf1*(rhof2 - rhof1)*(g_ & mesh_.Sf())
        )
    );

    //- Predicted fluxes
    surfaceScalarField phiHbyA1
    (
        IOobject::groupName("phiHbyA", fluid1_.name()),
        rAUf1
       *(
                alphaRhof10*MRF_.absolute(phi1.oldTime())/mesh_.time().deltaT()
            +   fvc::flux(UEqn1.H())
          
        )
        -   phig1
    );
    surfaceScalarField phiHbyA2
    (
        IOobject::groupName("phiHbyA", fluid2_.name()),
        rAUf2
       *(
                alphaRhof20*MRF_.absolute(phi2.oldTime())/mesh_.time().deltaT()
            +   fvc::flux(UEqn2.H())
          
        )
        -   phig2
    );

    if (!partialElimination_)
    {
        phiHbyA1 += rAUf1*Kdf12*MRF_.absolute(phi2);
        phiHbyA2 += rAUf2*Kdf12*MRF_.absolute(phi1);
    }

    //- Predicted mixture flux
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        alphaf1*phiHbyA1
      + alphaf2*phiHbyA2
    );

    if (partialElimination_)
    {
        phiHbyA += 
            alphaf1*rAUf1*Kdf12*MRF_.absolute(phi2) +
            alphaf2*rAUf2*Kdf12*MRF_.absolute(phi1);
    } 

    MRF_.makeRelative(phiHbyA);

    // Interpolate mixture 1/A field
    surfaceScalarField rAUf
    (
        "rAUf",
        mag(alphaf1*alpharAUf1 + alphaf2*alpharAUf2)
    );

    // Update the fixedFluxPressure BCs to ensure flux consistency
    surfaceScalarField::Boundary phib(phi_.boundaryField());
    phib = 
            alphaf1.boundaryField()*phi1.boundaryField() 
        +   alphaf2.boundaryField()*phi2.boundaryField();
    setSnGrad<fixedFluxPressureFvPatchScalarField>
    (
        p_rgh_.boundaryFieldRef(),
        (
            phiHbyA.boundaryField() - phib
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
                    ((phi2 + rAUf2*Kdf12*phi1) - (phi1 + rAUf1*Kdf12*phi2))
                   /(1 - rAUf1*rAUf2*sqr(Kdf12))
                );

                phi1 = phi_/movingAlphaf - alphaf2*phir;
                phi2 = phi_/movingAlphaf + alphaf1*phir;
            }

            fluid1_.dgdt() = -pEqnComp1 & p_rgh_;
            fluid2_.dgdt() = -pEqnComp2 & p_rgh_;
            
            p_rgh_.relax();

            mSfGradp = pEqnIncomp.flux()/rAUf;

            U1 = fvc::reconstruct(MRF_.absolute(phi1*movingAlphaf))/movingAlpha_;
            U2 = fvc::reconstruct(MRF_.absolute(phi2*movingAlphaf))/movingAlpha_;
            //U1 = fvc::reconstruct(MRF_.absolute(phi1));
            //U2 = fvc::reconstruct(MRF_.absolute(phi2));
            U1.correctBoundaryConditions();
            fvOptions_.correct(U1);            
            U2.correctBoundaryConditions();
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
