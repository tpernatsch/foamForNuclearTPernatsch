surfaceScalarField& phi1(fluid1_.phi());
surfaceScalarField& phi2(fluid2_.phi());
dimensionedScalar residualAlpha1(fluid1_.residualAlpha());
dimensionedScalar residualAlpha2(fluid2_.residualAlpha());
surfaceScalarField alphaf1(fvc::interpolate(alpha1));
surfaceScalarField alphaf2(fvc::interpolate(alpha2));
surfaceScalarField movingAlphaf(fvc::interpolate(movingAlpha_));

//- The faceMomentum algorithm relies on manually adding the time derivative
//  to the pressure equation, so that it can be defined on the faces (i.e.
//  deal with it in terms of fluxes). Thus, the time derivative coefficients
//  are defined here
surfaceScalarField alphaRhof1byDt
(
    "alphaRhof1",
    fvc::interpolate(max(alpha1, residualAlpha1)*rho1)/
    mesh_.time().deltaT()
);
surfaceScalarField alphaRhof2byDt
(
    "alphaRhof2",
    fvc::interpolate(max(alpha2, residualAlpha2)*rho2)/
    mesh_.time().deltaT()
);
surfaceScalarField alphaRhof10byDt
(
    "alphaRhof10",
    fvc::interpolate(max(alpha1.oldTime(), residualAlpha1)*rho1.oldTime())/
    mesh_.time().deltaT()
);
surfaceScalarField alphaRhof20byDt
(
    "alphaRhof20",
    fvc::interpolate(max(alpha2.oldTime(), residualAlpha2)*rho2.oldTime())/
    mesh_.time().deltaT()
);

surfaceScalarField Kdf12(fvc::interpolate(Kd12));

if(!experimentalPartialElimination_)
{
    //- Construct 1/A fields and limit to avoid problems for alpha = 0
    surfaceScalarField rAUf1
    (
        "rAUf."+fluid1_.name(),
        1.0/
        (
            fvc::interpolate(UEqn1.A())
        +   alphaRhof1byDt
        )
    );
    surfaceScalarField rAUf2
    (
        "rAUf."+fluid2_.name(),
        1.0/
        (
            fvc::interpolate(UEqn2.A())
        +   alphaRhof2byDt
        )
    );

    //- Interpolate 1/A to faces and limit
    surfaceScalarField alpharAUf1
    (
        max(alphaf1, residualAlpha1)*rAUf1
    );
    surfaceScalarField alpharAUf2
    (
        max(alphaf2, residualAlpha2)*rAUf2
    );

    while (pimple_.correct())
    {
        this->correctContErrs();

        rho_ = this->rho();

        //- Correct p_rgh for consistency with p and the updated densities
        p_rgh_ = p_ - rho_*gh_;

        surfaceScalarField deltaRhof12
        (
            fvc::interpolate(rho1)-fvc::interpolate(rho2)
        );

        //- Correct fixed-flux BCs to be consistent with the velocity BCs
        MRF_.correctBoundaryFlux(U1, phi1);
        MRF_.correctBoundaryFlux(U2, phi2);

        //- Compute buoyancy fluxes
        surfaceScalarField ghSnGradRho
        (
            "ghSnGradRho",
            ghf_*fvc::snGrad(rho_)*mesh_.magSf()
        );
        surfaceScalarField phig1
        (
            "phig."+fluid1_.name(),
            alpharAUf1*
            (
                ghSnGradRho
            -   alphaf2*deltaRhof12*(g_ & mesh_.Sf())
            )
        );
        surfaceScalarField phig2
        (
            "phig."+fluid2_.name(),
            alpharAUf2*
            (
                ghSnGradRho
            +   alphaf1*deltaRhof12*(g_ & mesh_.Sf())
            )
        );

        //- Predicted fluxes
        surfaceScalarField phiHbyA1
        (
            IOobject::groupName("phiHbyA", fluid1_.name()),
            rAUf1*
            (
                alphaRhof10byDt*MRF_.absolute(phi1.oldTime())
            +   fvc::flux(UEqn1.H())
            )
            -   phig1
        );
        surfaceScalarField phiHbyA2
        (
            IOobject::groupName("phiHbyA", fluid2_.name()),
            rAUf2*
            (
                alphaRhof20byDt*MRF_.absolute(phi2.oldTime())
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
        +   alphaf2*phiHbyA2
        );

        if (partialElimination_)
        {
            phiHbyA += 
                alphaf1*rAUf1*Kdf12*MRF_.absolute(phi2)
            +   alphaf2*rAUf2*Kdf12*MRF_.absolute(phi1);
        } 

        MRF_.makeRelative(phiHbyA);

        //- Interpolate mixture 1/A field
        surfaceScalarField rAUf
        (
            "rAUf",
            mag(alphaf1*alpharAUf1 + alphaf2*alpharAUf2)
        );

        //- Update the fixedFluxPressure BCs to ensure flux consistency
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

                if (partialElimination_)
                {
                    surfaceScalarField phir
                    (
                        (
                            (phi2 + rAUf2*Kdf12*phi1) 
                        -   (phi1 + rAUf1*Kdf12*phi2)
                        )/(1 - rAUf1*rAUf2*sqr(Kdf12))
                    );
                    
                    phi1 = (phi_ - alphaf2*phir)/movingAlphaf;
                    phi2 = (phi_ + alphaf1*phir)/movingAlphaf;
                }

                fluid1_.dgdt() = -pEqnComp1 & p_rgh_;
                fluid2_.dgdt() = -pEqnComp2 & p_rgh_;
                
                p_rgh_.relax();

                mSfGradp = pEqnIncomp.flux()/rAUf;

                //- The ifs are just to avoid useless field reconstructions
                //  for integer values of the oscillationLimiterFraction
                if (oscillationLimiterFraction_ != 0.0)
                {
                    if (oscillationLimiterFraction_ == 1.0)
                    {
                        U1 = 
                            fvc::reconstruct(MRF_.absolute(phi1*movingAlphaf))/
                            movingAlpha_;
                        U2 = 
                            fvc::reconstruct(MRF_.absolute(phi2*movingAlphaf))/
                            movingAlpha_;
                    }
                    else
                    {
                        U1 = 
                            oscillationLimiterFraction_*
                            (
                                fvc::reconstruct
                                (
                                    MRF_.absolute(phi1*movingAlphaf)
                                )/
                                movingAlpha_
                            )
                        +   (1-oscillationLimiterFraction_)*
                            (   
                                fvc::reconstruct(MRF_.absolute(phi1))
                            );
                        U2 = 
                            oscillationLimiterFraction_*
                            (
                                fvc::reconstruct
                                (
                                    MRF_.absolute(phi2*movingAlphaf)
                                )/
                                movingAlpha_
                            )
                        +   (1-oscillationLimiterFraction_)*
                            (   
                                fvc::reconstruct(MRF_.absolute(phi2))
                            );
                        }
                }
                else
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
else //- i.e. if doing experimentalPartialElimination
{
    //- Diagonal coefficients with time derivative contribution
    surfaceScalarField Af1
    (
        fvc::interpolate(UEqn1.A())
    +   alphaRhof1byDt
    );
    surfaceScalarField Af2
    (
        fvc::interpolate(UEqn2.A())
    +   alphaRhof2byDt
    );

    //- Construct inverse coefficients. The additional terms come from
    //  partial elimination
    surfaceScalarField rAUf1
    (
        "rAUf."+fluid1_.name(),
        1.0/(Af1 - sqr(Kdf12)/Af2)
    );
    surfaceScalarField rAUf2
    (
        "rAUf."+fluid2_.name(),
        1.0/(Af2 - sqr(Kdf12)/Af1)
    );

    //-
    surfaceScalarField KdbyAf1
    (
        Kdf12/Af1
    );
    surfaceScalarField KdbyAf2
    (
        Kdf12/Af2
    );

    //- 
    surfaceScalarField limitedAlphaf1(max(alphaf1, residualAlpha1));
    surfaceScalarField limitedAlphaf2(max(alphaf2, residualAlpha2));
    surfaceScalarField alpharAUf1
    (
        rAUf1*(limitedAlphaf1 + limitedAlphaf2*KdbyAf2)
    );
    surfaceScalarField alpharAUf2
    (
        rAUf2*(limitedAlphaf2 + limitedAlphaf1*KdbyAf1)
    );

    while (pimple_.correct())
    {
        this->correctContErrs();

        rho_ = this->rho();

        //- Correct p_rgh for consistency with p and the updated densities
        p_rgh_ = p_ - rho_*gh_;

        surfaceScalarField deltaRhof12
        (
            fvc::interpolate(rho1) - fvc::interpolate(rho2)
        );

        //- Correct fixed-flux BCs to be consistent with the velocity BCs
        MRF_.correctBoundaryFlux(U1, phi1);
        MRF_.correctBoundaryFlux(U2, phi2);

        //- Add explicit contributions (i.e. gravity and explicit time 
        //  derivative)
        surfaceScalarField ghSnGradRho
        (
            "ghSnGradRho",
            ghf_*fvc::snGrad(rho_)*mesh_.magSf()
        );
        surfaceScalarField source1
        (
            //- Gravity/buoyancy
            -   alphaf1*
                (
                    ghSnGradRho
                -   alphaf2/movingAlphaf*
                    deltaRhof12*(g_ & mesh_.Sf())
                )
            //- Explicit part of time derivative   
            +   alphaRhof10byDt*MRF_.absolute(phi1.oldTime())
        );
        surfaceScalarField source2
        (
            //- Gravity/buoyancy
            -   alphaf2*
                (
                    ghSnGradRho
                -   alphaf1/movingAlphaf*
                -   deltaRhof12*(g_ & mesh_.Sf())
                )
            //- Explicit part of time derivative
            +   alphaRhof20byDt*MRF_.absolute(phi2.oldTime())
        );
        surfaceScalarField phiSource1
        (
            rAUf1*
            (
                source1
            +   KdbyAf2*source2
            )
        ); 
        surfaceScalarField phiSource2
        (
            rAUf2*
            (
                source2
            +   KdbyAf1*source1
            )
        );

        //- Predicted fluxes
        surfaceScalarField phiH1(fvc::flux(UEqn1.H()));
        surfaceScalarField phiH2(fvc::flux(UEqn2.H()));
        surfaceScalarField phiHbyA1
        (
            IOobject::groupName("phiHbyA", fluid1_.name()),
            rAUf1*(phiH1 + phiH2*KdbyAf2)
        +   phiSource1
        );
        surfaceScalarField phiHbyA2
        (
            IOobject::groupName("phiHbyA", fluid2_.name()),
            rAUf2*(phiH2 + phiH1*KdbyAf1)
        +   phiSource2
        );

        //- Predicted mixture flux
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

        //- Construct the compressibility parts of the pressure equation
        tmp<fvScalarMatrix> p_rghCorr(correction(fvm::ddt(p_rgh_)));
        tmp<fvScalarMatrix> pEqnComp1;
        tmp<fvScalarMatrix> pEqnComp2;
        pEqnComp1 =
            (
                fluid1_.contErr() // +dmdt in contErr
              - fvc::Sp(fvc::ddt(alpha1) + fvc::div(fluid1_.alphaPhi()), rho1)
            )/rho1
        +   (alpha1*fluid1_.thermo().psi()/rho1)*p_rghCorr();

        pEqnComp2 =
            (
                fluid2_.contErr() // -dmdt in contErr
              - fvc::Sp(fvc::ddt(alpha2) + fvc::div(fluid2_.alphaPhi()), rho2)
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
                
                p_rgh_.relax();

                //- The ifs are just to avoid useless field reconstructions
                //  for integer values of the oscillationLimiterFraction
                if (oscillationLimiterFraction_ != 0.0)
                {
                    if (oscillationLimiterFraction_ == 1.0)
                    {
                        U1 = 
                            fvc::reconstruct(MRF_.absolute(phi1*movingAlphaf))/
                            movingAlpha_;
                        U2 = 
                            fvc::reconstruct(MRF_.absolute(phi2*movingAlphaf))/
                            movingAlpha_;
                    }
                    else
                    {
                        U1 = 
                            oscillationLimiterFraction_*
                            (
                                fvc::reconstruct
                                (
                                    MRF_.absolute(phi1*movingAlphaf)
                                )/movingAlpha_
                            )
                        +   (1-oscillationLimiterFraction_)*
                            (   
                                fvc::reconstruct(MRF_.absolute(phi1))
                            );
                        U2 = 
                            oscillationLimiterFraction_*
                            (
                                fvc::reconstruct
                                (
                                    MRF_.absolute(phi2*movingAlphaf)
                                )/movingAlpha_
                            )
                        +   (1-oscillationLimiterFraction_)*
                            (   
                                fvc::reconstruct(MRF_.absolute(phi2))
                            );
                        }
                }
                else
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
