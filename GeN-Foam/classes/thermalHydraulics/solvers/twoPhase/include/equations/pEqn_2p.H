surfaceScalarField& phi1(fluid1_.phi());
surfaceScalarField& phi2(fluid2_.phi());
dimensionedScalar residualAlpha1(fluid1_.residualAlpha());
dimensionedScalar residualAlpha2(fluid2_.residualAlpha());
surfaceScalarField alphaf1(fvc::interpolate(alpha1));
surfaceScalarField alphaf2(fvc::interpolate(alpha2));
surfaceScalarField movingAlphaf(fvc::interpolate(movingAlpha_));

bool noPartialElimination
(
    partialEliminationMode_ == partialEliminationMode::none
);
bool explicitPartialElimination
(
    partialEliminationMode_ == partialEliminationMode::legacy
);
bool implicitPartialElimination
(
    partialEliminationMode_ == partialEliminationMode::implicit
);

if (!implicitPartialElimination)
{
    //- Construct 1/A fields and limit to avoid problems for alpha = 0
    volScalarField rAU1
    (
        "rAU."+fluid1_.name(),
        1.0/
        (
            UEqn1.A()
        +   max(residualAlpha1 - alpha1, 0.0)*
            rho1/mesh_.time().deltaT()
        )
    );
    volScalarField rAU2
    (
        "rAU."+fluid2_.name(),
        1.0/
        (
            UEqn2.A()
        +   max(residualAlpha2 - alpha2, 0.0)*
            rho2/mesh_.time().deltaT()
        )
    );

    //- Interpolate 1/A to faces and limit
    surfaceScalarField alpharAUf1
    (
        fvc::interpolate(alpha1*rAU1)
    );
    surfaceScalarField alpharAUf2
    (
        fvc::interpolate(alpha2*rAU2)
    );

    // --- Pressure corrector loop
    while (pimple_.correct())
    {
        this->correctContErrs();
        rho_ = this->rho();
        p_rgh_ = p_ - rho_*gh_;

        //- Correct fixed-flux BCs to be consistent with the velocity BCs
        MRF_.correctBoundaryFlux(U1, phi1);
        MRF_.correctBoundaryFlux(U2, phi2);

        //- Construct H/A fields and limit to avoid problems for alpha = 0
        volVectorField HbyA1
        (
            "HbyA."+fluid1_.name(),
            rAU1*
            (
                UEqn1.H()
            +   max(residualAlpha1 - alpha1, 0.0)*
                rho1*U1.oldTime()/mesh_.time().deltaT()
            )
        );
            
        volVectorField HbyA2
        (
            "HbyA."+fluid2_.name(),
            rAU2*
            (
                UEqn2.H()
            +   max(residualAlpha2 - alpha2, 0.0)*
                rho2*U2.oldTime()/mesh_.time().deltaT()
            )
        );

        //- Construct fluxes due to explicit terms namely buoyancy/gravity 
        //  and the momentum source, if present
        surfaceScalarField ghSnGradRho
        (
            "ghSnGradRho",
            ghf_*fvc::snGrad(rho_)*mesh_.magSf()
        );
        surfaceScalarField deltaRhof12(fvc::interpolate(rho1 - rho2));
        surfaceScalarField phis1
        (
            alpharAUf1*
            (
                ghSnGradRho
            -   alphaf2*deltaRhof12*(g_ & mesh_.Sf())
            )
        );
        surfaceScalarField phis2
        (
            alpharAUf2*
            (
                ghSnGradRho
            +   alphaf1*deltaRhof12*(g_ & mesh_.Sf())
            )
        );
        if (structure_.momentumSourcePtr().valid())
        {
            surfaceScalarField fluxMS(fvc::flux(structure_.momentumSource()));
            phis1 -= alpharAUf1*fluxMS;
            phis2 -= alpharAUf2*fluxMS;
        }

        //- Construct phase fluxes
        surfaceScalarField phiHbyA1
        (
            "phiHbyA."+fluid1_.name(),
            fvc::flux(HbyA1)
        -   phis1
        );
        surfaceScalarField phiHbyA2
        (
            "phiHbyA."+fluid2_.name(),
            fvc::flux(HbyA2)
        -   phis2
        );

        //- Kd12/A1, Kd12/A2
        volScalarField rAUKd1 = rAU1*Kd12;
        volScalarField rAUKd2 = rAU2*Kd12;
        surfaceScalarField rAUKdf1(fvc::interpolate(rAUKd1));
        surfaceScalarField rAUKdf2(fvc::interpolate(rAUKd2));

        //- Add back explicit interphase drag contributions if not doing
        //  partial elimination
        if (noPartialElimination)
        {
            HbyA1 += rAUKd1*U2;
            HbyA2 += rAUKd2*U1;
            phiHbyA1 += rAUKdf1*MRF_.absolute(phi2);
            phiHbyA2 += rAUKdf2*MRF_.absolute(phi1);
        }
        
        //- Construct mixture flux
        surfaceScalarField phiHbyA
        (
            "phiHbyA",
            alphaf1*phiHbyA1
        +   alphaf2*phiHbyA2
        );
        if (explicitPartialElimination)
        {
            phiHbyA += 
                alphaf1*rAUKdf1*MRF_.absolute(phi2)
            +   alphaf2*rAUKdf2*MRF_.absolute(phi1);
        }

        //- Interpolate mixture 1/A field
        surfaceScalarField rAUf
        (
            "rAUf",
            mag(alphaf1*alpharAUf1 + alphaf2*alpharAUf2)
        );

        //- Update the fixedFluxPressure BCs to ensure flux consistency
        setSnGrad<fixedFluxPressureFvPatchScalarField>
        (
            p_rgh_.boundaryFieldRef(),
            (
                phiHbyA.boundaryField()
            -   (
                    alphaf1.boundaryField()*phi1.boundaryField()
                +   alphaf2.boundaryField()*phi2.boundaryField()
                )
            )/(mesh_.magSf().boundaryField()*rAUf.boundaryField())
        );

        //- Construct the compressibility parts of the pressure equation
        tmp<fvScalarMatrix> p_rghCorr(correction(fvm::ddt(p_rgh_)));
        tmp<fvScalarMatrix> pEqnComp1;
        tmp<fvScalarMatrix> pEqnComp2;
        pEqnComp1 =
            (
                fluid1_.contErr() // +dmdt in contErr
            -   fvc::Sp(fvc::ddt(alpha1) + fvc::div(fluid1_.alphaPhi()), rho1)
            )/rho1
        +   (alpha1*fluid1_.thermo().psi()/rho1)*p_rghCorr();

        pEqnComp2 =
            (
                fluid2_.contErr() // -dmdt in contErr
            -   fvc::Sp(fvc::ddt(alpha2) + fvc::div(fluid2_.alphaPhi()), rho2)
            )/rho2
        +   (alpha2*fluid2_.thermo().psi()/rho2)*p_rghCorr();

        //- Cache p prior to solve for density update
        volScalarField p_rgh_0(p_rgh_);

        //- Iterate over the pressure equation to correct for non-orthogonality
        while (pimple_.correctNonOrthogonal())
        {
            fvScalarMatrix pEqnIncomp
            (
                fvc::div(phiHbyA)
            -   fvm::laplacian(rAUf, p_rgh_)
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
                surfaceScalarField pEqnFlux(pEqnIncomp.flux());

                phi_ = phiHbyA + pEqnFlux;

                surfaceScalarField mSfGradp("mSfGradp", pEqnFlux/rAUf);
                
                phi1 = phiHbyA1 + alpharAUf1*mSfGradp;
                phi2 = phiHbyA2 + alpharAUf2*mSfGradp;

                if (explicitPartialElimination)
                {
                    surfaceScalarField phir
                    (
                        ((phi1 + rAUKdf1*phi2) - (phi2 + rAUKdf2*phi1))/
                        (1 - rAUKdf1*rAUKdf2)
                    );

                    phi1 = (phi_ + alphaf2*phir)/movingAlphaf;
                    phi2 = (phi_ - alphaf1*phir)/movingAlphaf;
                }

                fluid1_.dgdt() = -pEqnComp1 & p_rgh_;
                fluid2_.dgdt() = -pEqnComp2 & p_rgh_;

                if (momentumMode_ == momentumMode::cellCentered)
                {
                    p_rgh_.relax();

                    mSfGradp = pEqnIncomp.flux()/rAUf;

                    volVectorField Us1 = 
                        HbyA1 + fvc::reconstruct(alpharAUf1*mSfGradp - phis1);
                    volVectorField Us2 = 
                        HbyA2 + fvc::reconstruct(alpharAUf2*mSfGradp - phis2);

                    if (explicitPartialElimination)
                    {
                        volVectorField U
                        (
                            alpha1*(Us1 + rAUKd1*U2) + alpha2*(Us2 + rAUKd2*U1)
                        );

                        volVectorField Ur
                        (
                            ((1 - rAUKd2)*Us1 - (1 - rAUKd1)*Us2)/
                            (1 - rAUKd1*rAUKd2)
                        );

                        U1 = (U + alpha2*Ur)/movingAlpha_;
                        U2 = (U - alpha1*Ur)/movingAlpha_;
                    }
                    else
                    {
                        U1 = Us1;
                        U2 = Us2;
                    }
                }
                else    //- momentumMode_ == 
                        //  momentumMode::cellCenteredFaceReconstruction
                {
                    U1 = fvc::reconstruct(MRF_.absolute(phi1));
                    U2 = fvc::reconstruct(MRF_.absolute(phi2));
                }
                
                U1.correctBoundaryConditions();
                fvOptions_.correct(U1);
                fluid1_.magU() = mag(U1);

                U2.correctBoundaryConditions();
                fvOptions_.correct(U2);
                fluid2_.magU() = mag(U2);

                //- Update mixture velocity
                U_ = this->U();
                U_.correctBoundaryConditions();
            }
        }
        
        //- Update and limit the static pressure
        p_ = max(p_rgh_ + rho_*gh_, pMin_);

        //- Limit p_rgh
        p_rgh_ = p_ - rho_*gh_;

        //- Update densities from change in p_rgh
        rho1 += fluid1_.thermo().psi()*(p_rgh_ - p_rgh_0);
        rho2 += fluid2_.thermo().psi()*(p_rgh_ - p_rgh_0);

        //- Correct p_rgh for consistency with p and the updated densities
        rho_ = (rho1*alpha1 + rho2*alpha2)/movingAlpha_;
        p_rgh_ = p_ - rho_*gh_;
        p_rgh_.correctBoundaryConditions();
    }
}
else    //- i.e. if doing implicitPartialElimination
{
    volScalarField A1
    (
        UEqn1.A() 
    +   max(residualAlpha1 - alpha1, 0.0)*rho1/mesh_.time().deltaT()
    );
    volScalarField A2
    (
        UEqn2.A() 
    +   max(residualAlpha2 - alpha2, 0.0)*rho2/mesh_.time().deltaT()
    );
    volScalarField rAU1
    (
        "rAU."+fluid1_.name(),
        1.0/(A1-sqr(Kd12)/A2)
    );
    volScalarField rAU2
    (
        "rAU."+fluid2_.name(),
        1.0/(A2-sqr(Kd12)/A1)
    );

    //- 
    volScalarField limitedAlpha1(max(alpha1, residualAlpha1));
    volScalarField limitedAlpha2(max(alpha2, residualAlpha2));
    surfaceScalarField alpharAUf1
    (
        fvc::interpolate
        (
            rAU1*
            (
                limitedAlpha1
            +   limitedAlpha2*Kd12/A2
            )
        )
    );
    surfaceScalarField alpharAUf2
    (
        fvc::interpolate
        (
            rAU2*
            (
                limitedAlpha2
            +   limitedAlpha1*Kd12/A1
            )
        )
    );

    while (pimple_.correct())
    {
        this->correctContErrs();
        rho_ = this->rho();
        p_rgh_ = p_ - rho_*gh_;

        //- Correct fixed-flux BCs to be consistent with the velocity BCs
        MRF_.correctBoundaryFlux(U1, phi1);
        MRF_.correctBoundaryFlux(U2, phi2);

        //- Construct H/A fields and limit to avoid problems for alpha = 0
        tmp<volVectorField> H1(UEqn1.H());
        tmp<volVectorField> H2(UEqn2.H());
        volVectorField HbyA1
        (
            "HbyA."+fluid1_.name(),
            rAU1
           *(
                H1()
            +   H2()*Kd12/A2
            +   max(residualAlpha1 - alpha1, 0.0)*
                rho1*U1.oldTime()/mesh_.time().deltaT()
            )
        );
            
        volVectorField HbyA2
        (
            "HbyA."+fluid2_.name(),
            rAU2
           *(
                H2()
            +   H1()*Kd12/A1
            +   max(residualAlpha2 - alpha2, 0.0)*
                rho2*U2.oldTime()/mesh_.time().deltaT()
            )
        );
            
        //- Construct fluxes due to explicit source terms, namely 
        //  buoyancy/gravity effects and the momentum source, if it
        //  exists
        surfaceScalarField ghSnGradRho
        (
            "ghSnGradRho",
            ghf_*fvc::snGrad(rho_)*mesh_.magSf()
        );
        surfaceScalarField deltaRhof12(fvc::interpolate(rho1 - rho2));
        surfaceScalarField source1
        (
            alphaf1*
            (
                ghSnGradRho
            -   alphaf2/movingAlphaf*
                deltaRhof12*(g_ & mesh_.Sf())
            )
        );
        surfaceScalarField source2
        (
            alphaf2*
            (
                ghSnGradRho
            -   alphaf1/movingAlphaf*
            -   deltaRhof12*(g_ & mesh_.Sf())
            )
        );
        if (structure_.momentumSourcePtr().valid())
        {
            surfaceScalarField fluxMS(fvc::flux(structure_.momentumSource()));
            source1 -= alphaf1*fluxMS;
            source2 -= alphaf2*fluxMS;
        }
        surfaceScalarField phis1
        (
            fvc::interpolate(rAU1)*
            (
                source1
            +   fvc::interpolate(Kd12/A2)*source2
            )
        ); 
        surfaceScalarField phis2
        (
            fvc::interpolate(rAU2)*
            (
                source2
            +   fvc::interpolate(Kd12/A1)*source1
            )
        );

        //- Construct phase fluxes
        surfaceScalarField phiHbyA1
        (
            "phiHbyA."+fluid1_.name(),
            fvc::flux(HbyA1)
        -   phis1
        );
        surfaceScalarField phiHbyA2
        (
            "phiHbyA."+fluid2_.name(),
            fvc::flux(HbyA2)
        -   phis2
        );
        
        //- Construct the mean mixture flux
        surfaceScalarField phiHbyA
        (
            "phiHbyA",
            alphaf1*phiHbyA1
        +   alphaf2*phiHbyA2
        );

        MRF_.makeRelative(phiHbyA);

        //- Interpolate mixture 1/A field
        surfaceScalarField rAUf
        (
            "rAUf",
            mag(alphaf1*alpharAUf1 + alphaf2*alpharAUf2)
        );

        //- Update the fixedFluxPressure BCs to ensure flux consistency
        setSnGrad<fixedFluxPressureFvPatchScalarField>
        (
            p_rgh_.boundaryFieldRef(),
            (
                phiHbyA.boundaryField()
            -   (
                    alphaf1.boundaryField()*phi1.boundaryField()
                +   alphaf2.boundaryField()*phi2.boundaryField()
                )
            )/(mesh_.magSf().boundaryField()*rAUf.boundaryField())
        );

        //- Construct the compressibility parts of the pressure equation
        tmp<fvScalarMatrix> p_rghCorr(correction(fvm::ddt(p_rgh_)));
        tmp<fvScalarMatrix> pEqnComp1;
        tmp<fvScalarMatrix> pEqnComp2;
        
        pEqnComp1 =
            (
                fluid1_.contErr() // +dmdt in contErr
            -   fvc::Sp(fvc::ddt(alpha1) + fvc::div(fluid1_.alphaPhi()), rho1)
            )/rho1
        +   (alpha1*fluid1_.thermo().psi()/rho1)*p_rghCorr();

        pEqnComp2 =
            (
                fluid2_.contErr() // -dmdt in contErr
            -   fvc::Sp(fvc::ddt(alpha2) + fvc::div(fluid2_.alphaPhi()), rho2)
            )/rho2
        +   (alpha2*fluid2_.thermo().psi()/rho2)*p_rghCorr();

        //- Cache p prior to solve for density update
        volScalarField p_rgh_0(p_rgh_);

        //- Iterate over the pressure equation to correct for non-orthogonality
        while (pimple_.correctNonOrthogonal())
        {
            fvScalarMatrix pEqnIncomp
            (
                fvc::div(phiHbyA)
            -   fvm::laplacian(rAUf, p_rgh_)
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
                surfaceScalarField pEqnFlux(pEqnIncomp.flux());

                phi_ = phiHbyA + pEqnFlux;

                surfaceScalarField mSfGradp("mSfGradp", pEqnFlux/rAUf);
                
                phi1 = phiHbyA1 + alpharAUf1*mSfGradp;
                phi2 = phiHbyA2 + alpharAUf2*mSfGradp;

                fluid1_.dgdt() = -pEqnComp1 & p_rgh_;
                fluid2_.dgdt() = -pEqnComp2 & p_rgh_;

                if (momentumMode_ == momentumMode::cellCentered)
                {
                    p_rgh_.relax();

                    mSfGradp = pEqnIncomp.flux()/rAUf;

                    U1 = HbyA1 + fvc::reconstruct(alpharAUf1*mSfGradp - phis1);
                    U2 = HbyA2 + fvc::reconstruct(alpharAUf2*mSfGradp - phis2);
                }
                else    //- momentumMode_ == 
                        //  momentumMode::cellCenteredFaceReconstruction
                {
                    U1 = fvc::reconstruct(MRF_.absolute(phi1));
                    U2 = fvc::reconstruct(MRF_.absolute(phi2));
                }
                
                U1.correctBoundaryConditions();
                fvOptions_.correct(U1);
                fluid1_.magU() = mag(U1);
                
                U2.correctBoundaryConditions();
                fvOptions_.correct(U2);
                fluid2_.magU() = mag(U2);

                //- Update mixture velocity
                U_ = this->U();
                U_.correctBoundaryConditions();
            }
        }
        
        //- Update and limit the static pressure
        p_ = max(p_rgh_ + rho_*gh_, pMin_);

        //- Limit p_rgh
        p_rgh_ = p_ - rho_*gh_;

        //- Update densities from change in p_rgh
        rho1 += fluid1_.thermo().psi()*(p_rgh_ - p_rgh_0);
        rho2 += fluid2_.thermo().psi()*(p_rgh_ - p_rgh_0);
        
        //- Correct p_rgh for consistency with p and the updated densities
        rho_ = (rho1*alpha1 + rho2*alpha2)/movingAlpha_;
        p_rgh_ = p_ - rho_*gh_;
        p_rgh_.correctBoundaryConditions();
    }
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
