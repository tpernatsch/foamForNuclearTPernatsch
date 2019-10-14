{
    //- Field references
    volScalarField& he1(fluid1_.thermo().he());
    volScalarField& he2(fluid2_.thermo().he());
    volScalarField& T1(fluid1_.thermo().T());
    volScalarField& T2(fluid2_.thermo().T());
    const volVectorField& U1(fluid1_.U());
    const volVectorField& U2(fluid2_.U());
    const volScalarField& H1i(htcs_[fluid1_.name()+"."+fluid2_.name()]);
    const volScalarField& H2i(htcs_[fluid2_.name()+"."+fluid1_.name()]);
    const volScalarField& contErr1(fluid1_.contErr());
    const volScalarField& contErr2(fluid2_.contErr());
    volScalarField& rho1(fluid1_.thermo().rho());
    volScalarField& rho2(fluid2_.thermo().rho());
    tmp<volScalarField> cp1(fluid1_.thermo().Cp());
    tmp<volScalarField> cp2(fluid2_.thermo().Cp());

    //- Kinetic energy per unit volume
    volScalarField K1(0.5*U1&U1);
    volScalarField K2(0.5*U2&U2);

    //- Limiter fields
    volScalarField diffAlpha1(fluid1_-fluid1_.thermoResidualAlpha());
    volScalarField diffAlpha2(fluid2_-fluid2_.thermoResidualAlpha());
    volScalarField pos1(pos0(diffAlpha1));
    volScalarField pos2(pos0(diffAlpha2));
    volScalarField neg1(neg(diffAlpha1));
    volScalarField neg2(neg(diffAlpha2));

    //- Set interfacial temperatures and enthalpies/internal energies
    //  In particular, for no phaseChange, the interface is set
    //  so that heatFluxes across it are the same on both sides if
    //  both phases are present (pos1*pos2). Otherwise, it is set to
    //  Ti if phase i is present and phase j is not.
    if (phaseChange_)
    {
        scalar f(mesh_.fieldRelaxationFactor("iT12"));
        iT12_ = (1-f)*iT12_ + f*saturation_->Tsat(p_);
    }
    else
    {
        iT12_ = 
            pos1*pos2*(H1i*T1 + H2i*T2)
        /   (
                max
                (   
                    H1i+H2i, 
                    dimensionedScalar
                    (
                        "", 
                        dimPower/dimArea/dimTemperature, 
                        1e-6)
                    )
                )
        +   pos1*neg2*T1
        +   neg1*pos2*T2;
    }
    volScalarField he1i(fluid1_.thermo().he(p_, iT12_));
    volScalarField he2i(fluid2_.thermo().he(p_, iT12_));

    //- Avoid 0 matrix coeffs in the matrices, limit alphas
    volScalarField alpha1(fluid1_);
    volScalarField alpha2(fluid2_);
    alpha1 = max(alpha1, fluid1_.residualAlpha());
    alpha2 = max(alpha2, fluid2_.residualAlpha());

    volTensorField alphaEff1
    (
        structure_.tortuosity()
       *fluid1_.thermo().alphaEff(fluid1_.turbulence().mut())()
    );
    volTensorField alphaEff2
    (
        structure_.tortuosity()
       *fluid2_.thermo().alphaEff(fluid2_.turbulence().mut())()
    );

    //- Energy equations
    fvScalarMatrix E1Eqn
    (
            fvm::ddt(alpha1, rho1, he1) 
        +   fvm::div(fluid1_.alphaRhoPhi(), he1, "div(alphaRhoPhi,he)")
        -   fvm::Sp(contErr1, he1)
        +   fvc::ddt(alpha1, rho1, K1) 
        +   fvc::div(fluid1_.alphaRhoPhi(), K1, "div(alphaRhoPhi,K)")
        -   contErr1*K1
        +   alpha1*dpdt_
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
        -   fvm::Sp(contErr2, he2)
        +   fvc::ddt(alpha2, rho2, K2) 
        +   fvc::div(fluid2_.alphaRhoPhi(), K2, "div(alphaRhoPhi,K)")
        -   contErr2*K2
        +   alpha2*dpdt_
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

    if (phaseChange_)
    {
        //- Limit interfacial area so boiling can start 
        //  (very crude, good for now)
        iA12_ = max(iA12_, dimensionedScalar("", dimArea/dimVol, 1e-9));
        
        //- Compute the latent heat of vaporization in a manner consistent
        //  with how massTransfer related enthalpy is added to the energy
        //  equations. In particular, assuming 2 is vapour and 1 is liquid
        //  (but it works the other way around too):
        //  - when evaporating, mass is removed from the liquid at its current
        //    enthalpy and is added to the vapour at vapour saturation
        //    entalphy;
        //  - when condensing, mass is removed from the vapour at its current
        //    enthalpy and is added to the liquid at liquid saturation
        //    enthalpy.
        //    For example, in the case of condensing super-heated steam, its
        //    enthalpy will be larger than saturation enthalpy. However, if
        //    you remove mass from the steam at an enthalpy lower than its
        //    current enthalpy, the remaining steam will inevitably be
        //    hotter.
        //  The un-adjusted latent heat itself would be
        //  fluid2_.thermo().hc() - fluid1_.thermo().hc(), with hc() being
        //  the enthalpy of formation specified under the Hf keyword in the
        //  thermophysical dict of each phase
        volScalarField L
        (
                (
                        neg0(dmdt_)*(fluid2_.thermo().hc()) 
                    +   pos(dmdt_)*(he1 - he1i + fluid2_.thermo().hc())
                )
            -   (
                        pos0(dmdt_)*(fluid1_.thermo().hc()) 
                    +   neg(dmdt_)*(he2 - he2i + fluid1_.thermo().hc())
                )
        );

        scalar f(mesh_.fieldRelaxationFactor("dmdt"));
        if (pimple_.finalIter())
        {
            f = 1.0;
        }

        //- Compute mass transfer according to chosen model
        switch (phaseChangeModel_)
        {
            case phaseChangeModel::heatConductionLimited : 

                dmdt_ = 
                        (1.0-f)*dmdt_ 
                    +   f*((H1i*iA12_*(T1-iT12_)) + (H2i*iA12_*(T2-iT12_)))/L;
                break;

            case phaseChangeModel::individualSideDriven :
                dmdt_ = 
                        (1.0-f)*dmdt_ 
                    +   
                        f*
                        (
                            (fluid2_.name() == volatilePhaseName_) ?
                            posPart(H1i*iA12_*(T1-iT12_)/L) + 
                            negPart(H2i*iA12_*(T2-iT12_)/L)
                            :
                            negPart(H1i*iA12_*(T1-iT12_)/L) + 
                            posPart(H2i*iA12_*(T2-iT12_)/L)
                        );
                break;

            default : break;
        }
        
        //- Prevent removing mass from phases that are not present in a cell
        dmdt_ = posPart(dmdt_)*pos(fluid1_) + negPart(dmdt_)*pos(fluid2_);

        //- Add energy contributions as explained above
        volScalarField dmdt12(posPart(dmdt_));
        volScalarField dmdt21(-negPart(dmdt_));
        E1Eqn -= dmdt21*he1i - fvm::Sp(dmdt12, he1);
        E2Eqn -= dmdt12*he2i - fvm::Sp(dmdt21, he2);

        Info<< "L (avg min max) = "
            << " " << min(L).value()
            << " " << max(L).value()
            << " J/kg" << endl;
        Info<< "dmdt " << fluid1_.name() << "->" << fluid2_.name() 
            << " (avg min max) = " << dmdt_.weightedAverage(mesh_.V()).value()
            << " " << min(dmdt_).value()
            << " " << max(dmdt_).value()
            << " kg/m3/s" << endl;
    }
    
    //- Add linearized interfacial heat transfer contributions
    E1Eqn -= 
        (
            H1i*iA12_*(he1/cp1() + iT12_ - T1)
        -   fvm::Sp(H1i*iA12_/cp1(), he1)
        );
    E2Eqn -= 
        (
            H2i*iA12_*(he2/cp2() + iT12_ - T2)
        -   fvm::Sp(H2i*iA12_/cp2(), he2)
        );

    //- Structure heat sources
    //- Correct surface temperature of the power-generating subStructures 
    //  of the structures and add source/sink terms to the fluids energy 
    //  eqns
    volScalarField H1sfrac(htcs_[fluid1_.name()+".structure"]*frac1_);
    volScalarField H2sfrac(htcs_[fluid2_.name()+".structure"]*frac2_);
    structure_.correct
    (
        H1sfrac*T1 + H2sfrac*T2, 
        H1sfrac + H2sfrac
    );
    E1Eqn -= structure_.linearizedSemiImplicitHeatSource(fluid1_, H1sfrac)();
    E2Eqn -= structure_.linearizedSemiImplicitHeatSource(fluid2_, H2sfrac)();

    //- Relax and solve equations, re-limit enthalpies, update fluid
    //  temperatures
    E1Eqn.relax();
    fvOptions_.constrain(E1Eqn);
    E1Eqn.solve();
    fvOptions_.correct(he1);
    he1 = 
        max
        (
            pos1*he1 + neg1*he1i, 
            dimensionedScalar("", dimEnergy/dimMass, 100)
        );
    fluid1_.thermo().correct();

    E2Eqn.relax();
    fvOptions_.constrain(E2Eqn);
    E2Eqn.solve();
    fvOptions_.correct(he2);
    he2 = 
        max
        (
            pos2*he2 + neg2*he2i, 
            dimensionedScalar("", dimEnergy/dimMass, 100)
        );
    fluid2_.thermo().correct();

    //- Energy conservation
    volScalarField dmdt12(posPart(dmdt_));
    volScalarField dmdt21(-negPart(dmdt_));
    fluid1_.heErr() =
        pos1*
        (
                fvc::ddt(alpha1, rho1, he1)
            +   fvc::div(fluid1_.alphaRhoPhi(), he1, "div(alphaRhoPhi,T)")
            -   contErr1*he1
            -   fvc::laplacian
                (
                    fvc::interpolate(alpha1)
                   *fvc::interpolate(alphaEff1),
                    he1
                )
            +   fvc::ddt(alpha1, rho1, K1) 
            +   fvc::div(fluid1_.alphaRhoPhi(), K1, "div(alphaRhoPhi,K)")
            -   contErr1*K1
            +   alpha1*dpdt_
            -   H1i*iA12_*(iT12_-T1)
            -   structure_.explicitHeatSource(fluid1_, H1sfrac)
            -   alpha1*rho1*(U1&g_)
            -   (dmdt21*he1i - dmdt12*he1)
        )/(alpha1*rho1*he1);
    fluid1_.heErr().correctBoundaryConditions();
    fluid2_.heErr() =
        pos2*
        (
                fvc::ddt(alpha2, rho2, he2)
            +   fvc::div(fluid2_.alphaRhoPhi(), he2, "div(alphaRhoPhi,T)")
            -   contErr2*he2
            -   fvc::laplacian
                (
                    fvc::interpolate(alpha2)
                   *fvc::interpolate(alphaEff2),
                    he2
                )
            +   fvc::ddt(alpha2, rho2, K2) 
            +   fvc::div(fluid2_.alphaRhoPhi(), K2, "div(alphaRhoPhi,K)")
            -   contErr2*K2
            +   alpha2*dpdt_
            -   H2i*iA12_*(iT12_-T2)
            -   structure_.explicitHeatSource(fluid2_, H2sfrac)
            -   alpha2*rho2*(U2&g_)
            -   (dmdt12*he2i - dmdt21*he2)
        )/(alpha2*rho2*he2);
    fluid2_.heErr().correctBoundaryConditions();

    Info<< "T interface (avg min max) = "
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
    Info<< "Energy error " << fluid1_.name() << "(avg min max) = "
        << fluid1_.heErr().weightedAverage(mesh_.V()).value()
        << " " << min(fluid1_.heErr()).value()
        << " " << max(fluid1_.heErr()).value()
        << " 1/s" << endl;
    Info<< "Energy error " << fluid2_.name() << "(avg min max) = "
        << fluid2_.heErr().weightedAverage(mesh_.V()).value()
        << " " << min(fluid2_.heErr()).value()
        << " " << max(fluid2_.heErr()).value()
        << " 1/s" << endl;
}
