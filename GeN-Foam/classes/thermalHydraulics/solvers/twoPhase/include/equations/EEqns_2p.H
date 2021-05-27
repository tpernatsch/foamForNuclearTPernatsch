{
    //- Relax heat transfer coefficients in table if requested
    this->relaxFieldPtrTable<volScalarField>
    (
        word("htc"), 
        htcs_
    );

    //- Field references
    volScalarField& he1(fluid1_.thermo().he());
    volScalarField& he2(fluid2_.thermo().he());
    volScalarField& T1(fluid1_.thermo().T());
    volScalarField& T2(fluid2_.thermo().T());
    const volVectorField& U1(fluid1_.U());
    const volVectorField& U2(fluid2_.U());
    const volScalarField& H1i(*htcs_[fluid1_.name()+"."+fluid2_.name()]);
    const volScalarField& H2i(*htcs_[fluid2_.name()+"."+fluid1_.name()]);
    const volScalarField& contErr1(fluid1_.contErr());
    const volScalarField& contErr2(fluid2_.contErr());
    volScalarField& rho1(fluid1_.rho());
    volScalarField& rho2(fluid2_.rho());
    tmp<volScalarField> Cp1(fluid1_.thermo().Cp());
    tmp<volScalarField> Cp2(fluid2_.thermo().Cp());
    
    /*
    Info<< "H1i (avg min max) = "
        << H1i.weightedAverage(mesh_.V()).value()
        << " " << min(H1i).value()
        << " " << max(H1i).value()
        << " W/m2/K" << endl;
    Info<< "H2i (avg min max) = "
        << H2i.weightedAverage(mesh_.V()).value()
        << " " << min(H2i).value()
        << " " << max(H2i).value()
        << " W/m2/K" << endl;
    */

    //- Kinetic energy per unit volume
    volScalarField K1(0.5*U1&U1);
    volScalarField K2(0.5*U2&U2);

    //- Limiter fields
    volScalarField diffAlpha1
    (
        fluid1_.normalized()
    -   fluid1_.thermoResidualAlpha()
    );
    volScalarField diffAlpha2
    (
        fluid2_.normalized()
    -   fluid2_.thermoResidualAlpha()
    );
    volScalarField pos1(pos0(diffAlpha1));
    volScalarField pos2(pos0(diffAlpha2));
    volScalarField neg1(neg(diffAlpha1));
    volScalarField neg2(neg(diffAlpha2));
    volScalarField pos1pos2(pos1*pos2);

    //- Set interfacial temperatures and enthalpies/internal energies. In
    //  particular, for no phaseChange, the interface is set so that heatFluxes
    //  across it are the same on both sides if both phases are present
    //  (via pos1*pos2 marker field). Otherwise, it is set to Ti if phase i is
    //  present and phase j is not.
    if (phaseChange_.valid())
    {
        phaseChange_->correctInterfacialT();
    }
    else
    {
        iT12_ = 
            pos1pos2*(H1i*T1 + H2i*T2)/
            (
                max
                (   
                    H1i+H2i, 
                    dimensionedScalar
                    (
                        "", 
                        dimPower/dimArea/dimTemperature, 
                        1e-6
                    )
                )
            )
        +   pos1*neg2*T1
        +   neg1*pos2*T2;
    }
    
    volScalarField he1i(fluid1_.thermo().he(p_, iT12_));
    volScalarField he2i(fluid2_.thermo().he(p_, iT12_));

    //- Avoid 0 matrix coeffs by limiting the alphas that will figure in the
    //  equations
    volScalarField alpha1(fluid1_);
    volScalarField alpha2(fluid2_);
    alpha1 = max(alpha1, fluid1_.residualAlpha());
    alpha2 = max(alpha2, fluid2_.residualAlpha());

    //- Check if structure power off criteria are met and turn off power if so
    structure_.checkPowerOff();

    volTensorField alphaEff1
    (
        structure_.tortuosity()
       *fluid1_.thermo().alphaEff(fluid1_.turbulence().alphat())()
    );
    volTensorField alphaEff2
    (
        structure_.tortuosity()
       *fluid2_.thermo().alphaEff(fluid2_.turbulence().alphat())()
    );

    autoPtr<volScalarField> dpdtPtr;
    if (fluid1_.thermo().dpdt() or fluid2_.thermo().dpdt())
    {
        //- Pressure power term limited to 10 bar/s
        dimensionedScalar maxDpdt
        (
            "", 
            dimPressure/dimTime, 
            pimple_.dict().lookupOrDefault<scalar>("maxDpdt", 1e6)
        );
        dpdtPtr.reset
        (
            new volScalarField
            (
                max
                (
                    min
                    (
                        fvc::ddt(p_), maxDpdt
                    ),
                    -maxDpdt
                )
            )
        );
    }
    
    //- Construct enthalpy equations
    fvScalarMatrix E1Eqn
    (
            fvm::ddt(alpha1, rho1, he1) 
        +   fvm::div(fluid1_.alphaRhoPhi(), he1, "div(alphaRhoPhi,he)")
        -   fvm::laplacian
            (
                fvc::interpolate(alpha1)
               *fvc::interpolate(alphaEff1),
                he1
            )
        ==  
            alpha1*rho1*(U1&g_) 
        +   fvOptions_(alpha1, rho1, he1)
    );
    fvScalarMatrix E2Eqn
    (
            fvm::ddt(alpha2, rho2, he2)
        +   fvm::div(fluid2_.alphaRhoPhi(), he2, "div(alphaRhoPhi,he)")
        -   fvm::laplacian
            (
                fvc::interpolate(alpha2)
               *fvc::interpolate(alphaEff2),
                he2
            )
        ==  
            alpha2*rho2*(U2&g_) 
        +   fvOptions_(alpha2, rho2, he2)
    );

    if (fluid1_.thermo().dpdt())
        E1Eqn -= alpha1*dpdtPtr();
    if (fluid2_.thermo().dpdt())
        E2Eqn -= alpha2*dpdtPtr();

    scalar c
    (
        pimple_.dict().lookupOrDefault<scalar>
        (
            "continuityErrorScaleFactor", 
            1.0
        )
    );

    switch (heContErrCompensationMode_)
    {
        case contErrCompensationMode::none :
        {
            break;
        }
        case contErrCompensationMode::Su :
        {
            E1Eqn -= c*fluid1_.contErr()*he1;
            E2Eqn -= c*fluid2_.contErr()*he2;
            break;
        }
        case contErrCompensationMode::Sp :
        {
            E1Eqn -= fvm::Sp(c*fluid1_.contErr(), he1);
            E2Eqn -= fvm::Sp(c*fluid2_.contErr(), he2);
            break;
        }
        case contErrCompensationMode::SuSp :
        {
            E1Eqn += fvm::SuSp(-c*fluid1_.contErr(), he1);
            E2Eqn += fvm::SuSp(-c*fluid2_.contErr(), he2);
            break;
        }
    }

    //- Add kinetic energy contributions
    if (runTime_.controlDict().lookupOrDefault("includeKineticEnergy", true))
    {
        volScalarField K1(0.5*U1&U1);
        E1Eqn +=
        (
            fvc::ddt(alpha1, rho1, K1) 
        +   fvc::div(fluid1_.alphaRhoPhi(), K1, "div(alphaRhoPhi,K)")
        -   contErr1*K1
        );
        volScalarField K2(0.5*U2&U2);
        E2Eqn +=
        (
            fvc::ddt(alpha2, rho2, K2) 
        +   fvc::div(fluid2_.alphaRhoPhi(), K2, "div(alphaRhoPhi,K)")
        -   contErr2*K2
        );
    }

    //- Add mass transfer enthalpy contributions
    if (phaseChange_.valid())
    {
        //- Compute/update dmdt_
        phaseChange_->correctMassTransfer();

        //- Set phaseChange marker fields (currently only used by the
        //  nucleateBoiling model)
        scalarField& pCS1(fluid1_.phaseChangeSign());
        scalarField& pCS2(fluid2_.phaseChangeSign());
        forAll(dmdt_, i)
        {
            const scalar& dmdti(dmdt_[i]);
            pCS1[i] = (dmdti > 0) ? -1.0 : ((dmdti < 0) ? 1.0 : 0);
            pCS2[i] = -pCS1[i];
        }

        volScalarField dmdt12(posPart(dmdt_));
        volScalarField dmdt21(-negPart(dmdt_));

        E1Eqn -= dmdt21*he1i - fvm::Sp(dmdt12, he1);
        E2Eqn -= dmdt12*he2i - fvm::Sp(dmdt21, he2);

        Info<< "dmdt " << fluid1_.name() << "->" << fluid2_.name() 
            << " (avg min max) = " << dmdt_.weightedAverage(mesh_.V()).value()
            << " " << min(dmdt_).value()
            << " " << max(dmdt_).value()
            << " kg/m3/s" << endl;
    }
    
    //- Add linearized interfacial heat transfer contributions
    
    E1Eqn -= 
        (
            H1i*iA12_*(he1/Cp1() + iT12_ - T1)
        -   fvm::Sp(H1i*iA12_/Cp1(), he1)
        );
    E2Eqn -= 
        (
            H2i*iA12_*(he2/Cp2() + iT12_ - T2)
        -   fvm::Sp(H2i*iA12_/Cp2(), he2)
        );
    
    /*
    volScalarField noDmdt
    (
        pos0(dmdt_)*neg0(dmdt_)
    );
    E1Eqn -= 
        pos1pos2*
        (
            noDmdt*(
                H1i*iA12_*(he1/Cp1() + iT12_ - T1)
            -   fvm::Sp(H1i*iA12_/Cp1(), he1)
            )
        +   fvm::Sp(posPart(dmdt_)*phaseChange_->L()/he1, he1)
        );
    E2Eqn -= 
        pos1pos2*
        (
            noDmdt*(
                H2i*iA12_*(he2/Cp2() + iT12_ - T2)
            -   fvm::Sp(H2i*iA12_/Cp2(), he2)
            )
        +   pos(T2-iT12_)*fvm::Sp(negPart(dmdt_)*phaseChange_->L()/he2, he2)
        );
    */
    //- Structure heat sources
    //- Correct surface temperature of the power-generating subStructures 
    //  of the structures and add source/sink terms to the fluids energy 
    //  eqns. Also, limit fluid-structure HTCs only to regions where the fluid
    //  phase fraction is above thermoResidualAlpha
    volScalarField& H1s(*htcs_[fluid1_.name()+".structure"]);
    //H1s *= pos1;
    volScalarField& H2s(*htcs_[fluid2_.name()+".structure"]);
    //H2s *= pos2;
    
    E1Eqn -= structure_.linearizedSemiImplicitHeatSource(fluid1_, H1s)();
    E2Eqn -= structure_.linearizedSemiImplicitHeatSource(fluid2_, H2s)();

    //- Add additional sources deposited directly in the fluids, if any
    if (fluid1_.powerDensityPtr().valid())
    {
        E1Eqn -= fluid1_.powerDensity();
    }
    if (fluid2_.powerDensityPtr().valid())
    {
        E2Eqn -= fluid2_.powerDensity();
    }

    //- Add stabilization sources so that the temperatures will tend to their
    //  interfacial values if the phase fraction is below thermoResidualAlpha
    if (heStabilizationMode_ == heStabilizationMode::source)
    {
        dimensionedScalar HStab
        (
            "", 
            dimPower/dimVolume/dimTemperature, 
            1e10
        );
        if (fluid1_.thermoResidualAlpha().value() != 0.0)
        {
            volScalarField C1
            (
                (-neg1*diffAlpha1/fluid1_.thermoResidualAlpha())*HStab
            );
            E1Eqn -= 
            (
                C1*(he1/Cp1() + iT12_ - T1) 
            -   fvm::Sp(C1/Cp1(), he1)
            );
        }
        if (fluid2_.thermoResidualAlpha().value() != 0.0)
        {
            volScalarField C2
            (
                (-neg2*diffAlpha2/fluid2_.thermoResidualAlpha())*HStab
            );
            E2Eqn -= 
            (
                C2*(he2/Cp2() + iT12_ - T2) 
            -   fvm::Sp(C2/Cp2(), he2)
            );
        }
    }
    //- Store prev iter for structure T solution
    T1.storePrevIter();
    T2.storePrevIter();

    //- Relax and solve equations
    E1Eqn.relax();
    fvOptions_.constrain(E1Eqn);
    E1Eqn.solve();
    fvOptions_.correct(he1);
    he1.relax();
    E2Eqn.relax();
    fvOptions_.constrain(E2Eqn);
    E2Eqn.solve();
    fvOptions_.correct(he2);
    he2.relax(); 

    //- Stabilize enthalpies by cutoff if source mode not used
    if (heStabilizationMode_ == heStabilizationMode::cutoff)
    {
        he1 = 
            max
            (
                pos1*he1 + neg1*he1i, 
                dimensionedScalar("", dimEnergy/dimMass, 100)
            );
        he2 = 
            max
            (
                pos2*he2 + neg2*he2i, 
                dimensionedScalar("", dimEnergy/dimMass, 100)
            );
    }
    
    //- Correct temperatures based on new enthalpy and update other
    //  temperature-dpendent quantities
    fluid1_.thermo().correct();
    fluid2_.thermo().correct();

    //- Update structure temperature
    structure_.correct
    (
        H1s*T1.prevIter() + H2s*T2.prevIter(), 
        H1s + H2s
    );
    
    //- Correct fluid Prandtl numbers
    fluid1_.Pr() = Cp1()*fluid1_.thermo().mu()/fluid1_.thermo().kappa();
    fluid2_.Pr() = Cp2()*fluid2_.thermo().mu()/fluid2_.thermo().kappa();
    
    Info<< "T.interface (avg min max) = "
        << iT12_.weightedAverage(mesh_.V()).value()
        << " " << min(iT12_).value()
        << " " << max(iT12_).value()
        << " K" << endl;
    Info<< T1.name() << " (avg min max) = "
        << T1.weightedAverage(mesh_.V()).value()
        << " " << min(T1).value()
        << " " << max(T1).value()
        << " K" << endl;
    Info<< T2.name() << " (avg min max) = "
        << fluid2_.T().weightedAverage(mesh_.V()).value()
        << " " << min(T2).value()
        << " " << max(T2).value()
        << " K" << endl;
}
