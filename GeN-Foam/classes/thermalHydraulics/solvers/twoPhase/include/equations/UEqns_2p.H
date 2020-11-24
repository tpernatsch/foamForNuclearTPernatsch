//- Field refs
volVectorField& U1(fluid1_.U());
volVectorField& U2(fluid2_.U());
const volScalarField& alpha1(fluid1_);
const volScalarField& alpha2(fluid2_);
volScalarField& rho1(fluid1_.thermo().rho());
volScalarField& rho2(fluid2_.thermo().rho());
fvVectorMatrix UEqn1(U1, dimMass*dimLength/dimTime/dimTime);
fvVectorMatrix UEqn2(U2, dimMass*dimLength/dimTime/dimTime);
volTensorField& Kd12t(*Kds_[fluid1_.name()+"."+fluid2_.name()]);
volTensorField& Kd1st(*Kds_[fluid1_.name()+".structure"]);
volTensorField& Kd2st(*Kds_[fluid2_.name()+".structure"]);

//- Correct turbulence. What about this copy-setting fuss? Well, the standard
//  OpenFOAM turbulence models (e.g. kEpsilon) do not limit the phase fraction
//  before constructing the epsilon and k equations. This means that if a phase
//  fraction happens to be 0 somewhere, it will end up causing some null
//  entries in the k or epsilon matrices (pure genius, whomever made it...).
//  The solutions are: 1) copy some of the standard kEpslion models and make
//  them part of FFSEulerFoam and make the changes directly in the turbulence
//  classes; 2) limit the phase fractions before solving the turbulence models
//  and then reset them to the un-limited values cached before. Option 1)
//  is a NO-NO as I don't want to break compatibility with standard OpenFOAM
//  stuff and I want to avoid code duplication too (as if "code duplication" 
//  is a concern in a code that is 22k+ lines without even counting the 
//  underlying OpenFOAM library...)
volScalarField fluid10(fluid1_);
volScalarField fluid20(fluid2_);
fluid1_.volScalarField::operator=(max(fluid1_, fluid1_.residualAlpha()));
fluid2_.volScalarField::operator=(max(fluid2_, fluid2_.residualAlpha()));
fluid1_.turbulence().correct();
fluid2_.turbulence().correct();
fluid1_.volScalarField::operator=(fluid10);
fluid2_.volScalarField::operator=(fluid20);

//- Limit values of Kd
scalar KdFF0
(
    pimple_.dict().lookupOrDefault<scalar>("minKdFF", 1)
);
scalar KdFF1
(
    pimple_.dict().lookupOrDefault<scalar>("maxKdFF", KdFF0)
);
scalar a1
(
    pimple_.dict().lookupOrDefault<scalar>("maxAlphaKdFF", 0)
);
forAll(mesh_.C(), i)
{
    tensor& Kd12ti(Kd12t[i]);
    scalar Kdlim(KdFF0);
    if (KdFF1 != KdFF0 and a1 != 0)
    {
        const scalar& a(min(fluid1_.normalized()[i], fluid2_.normalized()[i]));
        if (a < a1)
            Kdlim = KdFF1;
    }
    Kd12ti[0] = max(Kd12ti[0], Kdlim);
    Kd12ti[4] = max(Kd12ti[4], Kdlim);
    Kd12ti[8] = max(Kd12ti[8], Kdlim);  
}
Kd12t.correctBoundaryConditions();
scalar KdFS
(
    pimple_.dict().lookupOrDefault<scalar>("minKdFS", 0)
);
if (KdFS != 0)
{
    forAll(mesh_.C(), i)
    {
        for (int j=0; j<3; j++)
        {
            scalar& Kd1sij(Kd1st[i][4*j]);
            Kd1sij = max(Kd1sij, KdFS);
            scalar& Kd2sij(Kd2st[i][4*j]);
            Kd2sij = max(Kd2sij, KdFS);
        }
    }
    Kd1st.correctBoundaryConditions();
    Kd2st.correctBoundaryConditions();
}

//- Relax drag coefficients in table if requested
this->relaxFieldPtrTable<volTensorField>
(
    word("Kd"), 
    Kds_
);

//- Drag between fluids is always a scalar, so I construct it from the 
//  first tensor component
volScalarField Kd12
(
    IOobject
    (
        "Kd12",
        mesh_.time().timeName(),
        mesh_
    ),
    Kd12t.component(0),
    zeroGradientFvPatchScalarField::typeName
);

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
    &   Kd1st
    &   structure_.Rg2l()
    ),
    zeroGradientFvPatchScalarField::typeName
);

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
    &   Kd2st
    &   structure_.Rg2l()
    ),
    zeroGradientFvPatchScalarField::typeName
);

/*
Info<< "Kd12 (avg min max) = "
        << Kd12.weightedAverage(mesh_.V()).value()
        << " " << min(Kd12).value()
        << " " << max(Kd12).value()
        << " kg/m3/s" << endl;
*/

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
            //- Implicit diagonal contribution of fluid-structure drag,
            //  explicit off-diagonal contribution
        +   fvm::Sp((1.0/3.0)*tr(Kd2s), U2) + (dev(Kd2s)&U2) 
        ==
            dmdt12*U1
        +   fvOptions_(alpha2, rho2, U2)
    );

    scalar c
    (
        pimple_.dict().lookupOrDefault<scalar>
        (
            "continuityErrorScaleFactor", 
            1.0
        )
    );

    switch (UContErrCompensationMode_)
    {
        case contErrCompensationMode::none :
        {
            break;
        }
        case contErrCompensationMode::Su :
        {
            UEqn1 -= c*fluid1_.contErr()*U1;
            UEqn2 -= c*fluid2_.contErr()*U2;
            break;
        }
        case contErrCompensationMode::Sp :
        {
            UEqn1 -= fvm::Sp(c*fluid1_.contErr(), U1);
            UEqn2 -= fvm::Sp(c*fluid2_.contErr(), U2);
            break;
        }
        case contErrCompensationMode::SuSp :
        {
            UEqn1 += fvm::SuSp(-c*fluid1_.contErr(), U1);
            UEqn2 += fvm::SuSp(-c*fluid2_.contErr(), U2);
            break;
        }
    }

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

    //- Vritual mass contribution, handled outside the regime framework for now
    if (Vm_ > 0)
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
            Vm_*
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
        scalar f(myOps::relaxationFactor(mesh_, "virtualMass"));
        if (f == 1.0) // i.e. if not under-relaxing
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
        else // i.e. if under-relaxing
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
