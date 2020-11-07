//- Field refs
volVectorField& U1(fluid1_.U());
volVectorField& U2(fluid2_.U());
const volScalarField& alpha1(fluid1_);
const volScalarField& alpha2(fluid2_);
volScalarField& rho1(fluid1_.thermo().rho());
volScalarField& rho2(fluid2_.thermo().rho());
fvVectorMatrix UEqn1(U1, dimMass*dimLength/dimTime/dimTime);
fvVectorMatrix UEqn2(U2, dimMass*dimLength/dimTime/dimTime);

//- Relax drag coefficients in table if requested
this->relaxFieldPtrTable<volTensorField>
(
    word("momentumTransfer"), 
    Kds_
);

//- Drag between fluids is always a scalar, so I convert it. Additionaly,
//  limit it to a minimum > 0 to avoid phase decoupling
volScalarField Kd12
(
    IOobject
    (
        "Kd12",
        mesh_.time().timeName(),
        mesh_
    ),
    (*Kds_[fluid1_.name()+"."+fluid2_.name()]).component(0),
    zeroGradientFvPatchScalarField::typeName
);
//- The Cd is the actual drag coefficient that is required only by some
//  turbulence models and it is never used elsewhere in the code. I compute it
//  here rather than in the turbulence models because, if I were to use Kd12,
//  the this reconstruction might not work as intended at small Kd12 values, as
//  Kd12 is limited to a minimum of residualKd (so the Cd calculated from
//  within the turbulence models might end up being severely over-estimated in
//  such scenarios)
Cd12_ = 
    (2.0)*Kd12*FFPair_.DhDispersed()/FFPair_.rhoContinuous()/
    max
    (
        fluid1_*fluid2_*FFPair_.magUr(), 
        dimensionedScalar("", dimVelocity, 1e-3)
    );

//- Limit to > 0 to avoid phase decoupling
Kd12 = 
    max
    (
        Kd12,
        dimensionedScalar("", dimDensity/dimTime, 1)
    );

//- Limit to a very high value in those regions that are single-phase only to
//  avoid numerical bullshit
volScalarField twoPhaseMarkerField
(
    pos(fluid1_-fluid1_.residualAlpha())*pos(fluid2_-fluid2_.residualAlpha())
);
Kd12 = 
    twoPhaseMarkerField*Kd12
+   (1.0-twoPhaseMarkerField)*dimensionedScalar("", dimDensity/dimTime, 1e5);
Kd12.correctBoundaryConditions(); //- <= REQUIRED OR SHIT HAPPENS IN MPI RUNS

/*
Info<< "Kd12 (avg min max) = "
        << Kd12.weightedAverage(mesh_.V()).value()
        << " " << min(Kd12).value()
        << " " << max(Kd12).value()
        << " kg/m3/s" << endl;
*/

//- Drag between fluid1 and structure. It is expressed in the local
//  structure frame coordinates, so rotate it back into the global
//  frame
volTensorField Kd1s
(
    IOobject
    (
        "Kd1s",
        mesh_.time().timeName(),
        mesh_
    ),
    (
        structure_.Rl2g() 
    &   (*Kds_[fluid1_.name()+".structure"])
    &   structure_.Rg2l()
    ),
    zeroGradientFvPatchScalarField::typeName
);
Kd1s.correctBoundaryConditions();

//- Drag between fluid2 and structure. It is expressed in the local
//  structure frame coordinates, so rotate it back into the global
//  frame
volTensorField Kd2s
(
    IOobject
    (
        "Kd2s",
        mesh_.time().timeName(),
        mesh_
    ),
    (
        structure_.Rl2g() 
    &   (*Kds_[fluid2_.name()+".structure"])
    &   structure_.Rg2l()
    ),
    zeroGradientFvPatchScalarField::typeName
);
Kd2s.correctBoundaryConditions();

//- Correct turbulence. What about this copy-setting fuss? Well, the standard
//  OpenFOAM turbulence models (e.g. kEpsilon) do not limit the phase fraction
//  before constructing the epsilon and k equations. This means that if a phase
//  fraction happens to be 0 somewhere, it will end up causing some null
//  entries in the k or epsilon matrices (pure genius, whomever made it...).
//  The solutions are: 1) copy some of the standard kEpslion models and make
//  them part of FFSEulerFoam and make the changes directly in the turbulence
//  classes; 2) limit the phase fractions before solving the turbulence models
//  and them re-setting them to the un-limited values cahced before. Option 1)
//  is a NO-NO as I don't want to break compatibility with standard OpenFOAM
//  stuff as well as avoiding code duplication... so, here's the damn ugly
//  solution 2)
volScalarField fluid10(fluid1_);
volScalarField fluid20(fluid2_);
fluid1_.volScalarField::operator=(max(fluid1_, fluid1_.residualAlpha()));
fluid2_.volScalarField::operator=(max(fluid2_, fluid2_.residualAlpha()));
fluid1_.turbulence().correct();
fluid2_.turbulence().correct();
fluid1_.volScalarField::operator=(fluid10);
fluid2_.volScalarField::operator=(fluid20);

{
    volTensorField nuEff1
    (
        structure_.tortuosity()
       *fluid1_.turbulence().nuEff()()
    );
    volTensorField nuEff2
    (
        structure_.tortuosity()
       *fluid2_.turbulence().nuEff()()
    );
    
    //- Mass transfer
    volScalarField dmdt12(posPart(dmdt_));
    volScalarField dmdt21(-negPart(dmdt_)); 
    
    UEqn1 =
    (
            fvm::div(fluid1_.alphaRhoPhi(), U1, "div(alphaRhoPhi,U)")
        -   fvm::SuSp(fluid1_.contErr(), U1)
        +   fvm::Sp(dmdt12, U1)
        //- The following is just a re-arrangement of div(nu*grad(U)),
        //  i.e. + fluid1_.turbulence().divDevRhoReff(U1). However, it is
        //  necessary to do this, otherwise the tensor nature of nuEff can't
        //  be resolved in anisotropic porous media
        +   (
                -   fvm::laplacian(alpha1*rho1*nuEff1, U1)
                -   fvc::div
                    (
                        alpha1*rho1*nuEff1&dev2(T(fvc::grad(U1))), 
                        "div(alphaRhoPhiNu,U)"
                    )
            )
            
        //+   fvm::Sp(Kd12, U1)
            //- Implicit diagonal contribution of fluid-structure drag, 
            //  explicit off-diagonal contribution
        +   fvm::Sp((1.0/3.0)*tr(Kd1s), U1) + (dev(Kd1s)&U1)
        ==
            dmdt21*U2
        +   fvOptions_(alpha1, rho1, U1)
    );
    UEqn2 =
    (
            fvm::div(fluid2_.alphaRhoPhi(), U2, "div(alphaRhoPhi,U)")
        -   fvm::SuSp(fluid2_.contErr(), U2)
        +   fvm::Sp(dmdt21, U2)
            //- The following is just a re-arrangement of div(nu*grad(U))
        +   (
                -   fvm::laplacian(alpha2*rho2*nuEff2, U2)
                -   fvc::div
                    (
                        alpha2*rho2*nuEff2&dev2(T(fvc::grad(U2))), 
                        "div(alphaRhoPhiNu,U)"
                    )
            )
        //+   fvm::Sp(Kd12, U2)
            //- Implicit diagonal contribution of fluid-structure drag,
            //  explicit off-diagonal contribution
        +   fvm::Sp((1.0/3.0)*tr(Kd2s), U2) + (dev(Kd2s)&U2) 
        ==
            dmdt12*U1
        +   fvOptions_(alpha2, rho2, U2)
    );

    //- The time derivative is manually added to the pressure equation
    //  contribution if doing faceMomentum. Otherwise, do it here
    if (momentumMode_ != momentumMode::faceCentered) 
    {
        UEqn1 += fvm::ddt(alpha1, rho1, U1);
        UEqn2 += fvm::ddt(alpha2, rho2, U2);
    }

    //
    if (partialEliminationMode_ != partialEliminationMode::implicit)
    {
        UEqn1 += fvm::Sp(Kd12, U1);
        UEqn2 += fvm::Sp(Kd12, U2);
    }

    //- Vritual mass coefficient, not handled by regimes for now
    dimensionedScalar Vm
    (
        this->lookupOrDefault<scalar>("virtualMassCoeff", 0.0)
    );
    if (Vm.value() > 0)
    {
        const surfaceScalarField& phi1(fluid1_.phi());
        const surfaceScalarField& phi2(fluid2_.phi());

        fvVectorMatrix VmVolForce1(U1, dimMass*dimLength/dimTime/dimTime);
        fvVectorMatrix VmVolForce2(U2, dimMass*dimLength/dimTime/dimTime);

        bool virtualMassUsesMixtureDensity
        (
            pimple_.dict().lookupOrDefault<bool>
            (
                "virtualMassUsesMixtureDensity", 
                false
            )
        );

        volScalarField alphaRhoVm
        (
            Vm*
            FFPair_.alphaDispersed()*
            (
                (virtualMassUsesMixtureDensity) ?
                (
                    fluid1_.normalized()*fluid1_.thermo().rho() 
                +   fluid2_.normalized()*fluid2_.thermo().rho()
                )() :
                FFPair_.rhoContinuous()
            )
        );
        volVectorField DDtU1 = 
            fvc::ddt(U1)
        +   fvc::div(phi1, U1)
        -   fvc::div(phi1)*U1;
        volVectorField DDtU2 =
            fvc::ddt(U2)
        +   fvc::div(phi2, U2)
        -   fvc::div(phi2)*U2;
        scalar f(1.0);
        if (mesh_.relaxField("virtualMass"))
            f = mesh_.fieldRelaxationFactor("virtualMass");
        if (f == 1.0)
        {
            UEqn1 +=  
                alphaRhoVm*
                (
                    fvm::ddt(U1)
                +   fvm::div(phi1, U1)
                +   fvm::SuSp(-fvc::div(phi1), U1)
                -   DDtU2
                );
            UEqn2 += 
                alphaRhoVm*
                (
                    fvm::ddt(U2)
                +   fvm::div(phi2, U2)
                +   fvm::SuSp(-fvc::div(phi2), U2)
                -   DDtU1
                );
        }
        else
        {
            VmVolForce1_ = 
                f*
                (
                    alphaRhoVm*
                    (
                        fvm::ddt(U1)
                    +   fvm::div(phi1, U1)
                    +   fvm::SuSp(-fvc::div(phi1), U1)
                    -   DDtU2
                    )
                )
            +   (1.0-f)*VmVolForce1_;
            VmVolForce2_ = 
                f*
                (
                    alphaRhoVm*
                    (
                        fvm::ddt(U2)
                    +   fvm::div(phi2, U2)
                    +   fvm::SuSp(-fvc::div(phi2), U2)
                    -   DDtU1
                    )
                )
            +   (1.0-f)*VmVolForce2_;
            UEqn1 += VmVolForce1_;
            UEqn2 += VmVolForce2_;
        }

        //- Modify MRF and add contributions
        MRF_.addAcceleration(U1, DDtU1);
        MRF_.addAcceleration(U2, DDtU2);
        UEqn1 += MRF_.DDt(alpha1*rho1 + alphaRhoVm, U1);
        UEqn2 += MRF_.DDt(alpha2*rho2 + alphaRhoVm, U2);
    }
    else //- Just add unmodified MRF contributions
    {
        UEqn1 += MRF_.DDt(alpha1*rho1, U1);
        UEqn2 += MRF_.DDt(alpha2*rho2, U2);
    }

    UEqn1.relax();
    fvOptions_.constrain(UEqn1);
    fvOptions_.correct(U1);
    UEqn2.relax();
    fvOptions_.constrain(UEqn2);
    fvOptions_.correct(U2);
}
