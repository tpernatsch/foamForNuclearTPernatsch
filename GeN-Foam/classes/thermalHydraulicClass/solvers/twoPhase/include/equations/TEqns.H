{
    /*
        NOTE: Why not solve for the enthalpy/internal energy like most OpenFOAM
        solvers do? Why am I solving directly for temperature?

        Mainly one reason: if I solve for temperature I can implicitly treat
        the phase-interface heat fluxes which is a considerable benefit, as
        these tend to be the absolute major source of instability in two-phase
        solvers, even more than mass-transfer-assosciated energy transfer.
        What is the drawback? Well, if the cp has a significant dependence 
        on temperature, I might not be as conservative as with the 
        traditional approach IF I am trying to resolve sharp temperature
        changes in time. No big deal as one can always reduce the time-step.
    */

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
    volScalarField cp1(fluid1_.thermo().Cp());
    volScalarField cp2(fluid2_.thermo().Cp());
    volScalarField TRef
    (
        IOobject
        (
            "TRef",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("", dimTemperature, mesh_.fieldRelaxationFactor("TRef")),
        zeroGradientFvPatchScalarField::typeName
    );
    he1 = fluid1_.thermo().he(p_, T1);
    he2 = fluid2_.thermo().he(p_, T2);

    //- Kinetic energy per unit volume
    volScalarField K1(0.5*U1&U1);
    volScalarField K2(0.5*U2&U2);

    //- Limiter fields
    //  These serve to set the temperatures to the interfacial
    //  temperature if the phaseFraction of a phase falls below
    //  thermoResidualAlpha
    dimensionedScalar thermoResidualAlpha(this->lookupOrDefault("thermoResidualAlpha", 1e-2));
    volScalarField diffAlpha1(fluid1_-thermoResidualAlpha);
    volScalarField diffAlpha2(fluid2_-thermoResidualAlpha);
    volScalarField pos1(pos0(diffAlpha1));
    volScalarField pos2(pos0(diffAlpha2));
    volScalarField neg1(neg(diffAlpha1));
    volScalarField neg2(neg(diffAlpha2));

    //- Limit alpha to avoid 0 matrix coefficients
    volScalarField alpha1(fluid1_);
    volScalarField alpha2(fluid2_);
    alpha1 = max(alpha1, fluid1_.residualAlpha());
    alpha2 = max(alpha2, fluid2_.residualAlpha());

    //- Interfacial temperature
    if (phaseChange_)
    {
        scalar f(mesh_.fieldRelaxationFactor("iT12"));
        iT12_ = (1-f)*iT12_ + f*saturation_->Tsat(p_);
        
        if (mesh_.time().timeIndex() - mesh_.time().startTimeIndex() == 1)
        {
            IOobject iT12Header
            (
                iT12_.name(),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ
            );
            if (!iT12_.typeHeaderOk<volScalarField>(true))
            {
                iT12_ = saturation_->Tsat(p_);
            }
        }
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
    volScalarField he1s(fluid1_.thermo().he(p_, iT12_));
    volScalarField he2s(fluid2_.thermo().he(p_, iT12_));

    fvScalarMatrix E1Eqn
    (
            cp1*fvm::ddt(alpha1, rho1, T1) 
        +   cp1*
            (
                fvm::div(fluid1_.alphaRhoPhi(), T1, "div(alphaRhoPhi,T)")
            -   fvm::SuSp(contErr1, T1)
            )
        -   fvm::laplacian
            (
                fvc::interpolate(alpha1)*fvc::interpolate
                (
                    fluid1_.thermo().kappaEff(fluid1_.turbulence().mut())
                ),
                T1
            )
        +   fvc::ddt(alpha1, rho1, K1) 
        +   fvc::div(fluid1_.alphaRhoPhi(), K1, "div(alphaRhoPhi,K)")
        -   contErr1*K1
        +   alpha1*dpdt_
        ==
            H1i*iA12_*(iT12_)
        -   fvm::Sp(H1i*iA12_, T1)
        +   alpha1*rho1*(U1&g_) 
        +   cp1*fvOptions_(alpha1, rho1, T1)
    );

    fvScalarMatrix E2Eqn
    (
            cp2*fvm::ddt(alpha2, rho2, T2)
        +   cp2*
            (
                fvm::div(fluid2_.alphaRhoPhi(), T2, "div(alphaRhoPhi,T)")
            -   fvm::SuSp(contErr2, T2)
            )

        -   fvm::laplacian
            (
                fvc::interpolate(alpha2)
               *fvc::interpolate
                (
                    fluid2_.thermo().kappaEff(fluid2_.turbulence().mut())
                ),
                T2
            )
        +   fvc::ddt(alpha2, rho2, K2) 
        +   fvc::div(fluid2_.alphaRhoPhi(), K2, "div(alphaRhoPhi,K)")
        -   contErr2*K2
        +   alpha2*dpdt_
        
        ==
            H2i*iA12_*(iT12_)
        -   fvm::Sp(H2i*iA12_, T2)
        +   alpha2*rho2*(U2&g_) 
        +   cp2*fvOptions_(alpha2, rho2, T2)
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
        volScalarField L
        (
                (
                        neg0(dmdt_)*(fluid2_.thermo().hc()) 
                    +   pos(dmdt_)*(he1 - he1s + fluid2_.thermo().hc())
                )
            -   (
                        pos0(dmdt_)*(fluid1_.thermo().hc()) 
                    +   neg(dmdt_)*(he2 - he2s + fluid1_.thermo().hc())
                )
        );
        //- Under-relaxed heat conduction limited model for mass transfer
        scalar f(mesh_.fieldRelaxationFactor("dmdt"));
        volScalarField sumQ
        (
            H1i*iA12_*(T1-iT12_) + H2i*iA12_*(T2-iT12_)
        );
        if (pimple_.finalIter())
        {
            f = 1.0;
        }
        dmdt_ = (1.0-f)*dmdt_ + f*sumQ/L;
        
        dmdt_ = posPart(dmdt_)*pos(fluid1_) + negPart(dmdt_)*pos(fluid2_);
        volScalarField dmdt12(posPart(dmdt_));
        volScalarField dmdt21(-negPart(dmdt_));

        //- Add energy contributions as explained above (to avoid run-aways)
        
        E1Eqn -= 
                dmdt21*(he1s-fluid1_.thermo().he(p_, iT12_*0)) 
            -   (dmdt12*(he1-cp1*T1) + fvm::Sp(dmdt12*cp1, T1));
        E2Eqn -= 
                dmdt12*(he2s-fluid2_.thermo().he(p_, iT12_*0)) 
            -   (dmdt21*(he2-cp2*T2) + fvm::Sp(dmdt21*cp2, T2));

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
    volScalarField Qs1(structure_.heatSource(T1, H1sfrac));
    volScalarField Qs2(structure_.heatSource(T2, H2sfrac));
    E1Eqn -= Qs1;
    E2Eqn -= Qs2;

    //- Relax and solve eqns
    E1Eqn.relax();
    fvOptions_.constrain(E1Eqn);
    E1Eqn.solve();
    fvOptions_.correct(T1);
    T1 = max(T1, dimensionedScalar("", dimTemperature, SMALL));
    T1 = pos1*T1 + neg1*iT12_;
    E2Eqn.relax();
    fvOptions_.constrain(E2Eqn);
    E2Eqn.solve();
    fvOptions_.correct(T2);
    T2 = max(T2, dimensionedScalar("", dimTemperature, SMALL));
    T2 = pos2*T2 + neg2*iT12_;

    he1 = fluid1_.thermo().he(p_, T1);
    he2 = fluid2_.thermo().he(p_, T2);

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
                   *fvc::interpolate
                    (
                        fluid1_.thermo().alphaEff(fluid1_.turbulence().mut())
                    ),
                    he1
                )
            +   fvc::ddt(alpha1, rho1, K1) 
            +   fvc::div(fluid1_.alphaRhoPhi(), K1, "div(alphaRhoPhi,K)")
            -   contErr1*K1
            +   alpha1*dpdt_
            +   H1i*iA12_*(T1-iT12_)
            -   Qs1
            -   alpha1*rho1*(U1&g_)
            -   (dmdt21*he1s - dmdt12*he1)
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
                   *fvc::interpolate
                    (
                        fluid2_.thermo().alphaEff(fluid2_.turbulence().mut())
                    ),
                    he2
                )
            +   fvc::ddt(alpha2, rho2, K2) 
            +   fvc::div(fluid2_.alphaRhoPhi(), K2, "div(alphaRhoPhi,K)")
            -   contErr2*K2
            +   alpha2*dpdt_
            +   H2i*iA12_*(T2-iT12_)
            -   Qs2
            -   alpha2*rho2*(U2&g_)
            -   (dmdt12*he2s - dmdt21*he2)
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
