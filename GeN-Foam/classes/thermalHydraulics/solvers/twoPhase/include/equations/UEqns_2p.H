//- Field refs
volVectorField& U1(fluid1_.U());
volVectorField& U2(fluid2_.U());
const volScalarField& alpha1(fluid1_);
const volScalarField& alpha2(fluid2_);
volScalarField& rho1(fluid1_.thermo().rho());
volScalarField& rho2(fluid2_.thermo().rho());
fvVectorMatrix UEqn1(U1, dimMass*dimLength/dimTime/dimTime);
fvVectorMatrix UEqn2(U2, dimMass*dimLength/dimTime/dimTime);

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
    (1.0/3.0)*tr(Kds_[fluid1_.name()+"."+fluid2_.name()]),
    zeroGradientFvPatchScalarField::typeName
);
Kd12 = max(Kd12, residualKd_);
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
        structure_.Rl2g() 
    &   Kds_[fluid1_.name()+".structure"]
    &   structure_.Rg2l()
);

//- Drag between fluid2 and structure. It is expressed in the local
//  structure frame coordinates, so rotate it back into the global
//  frame
volTensorField Kd2s
(
        structure_.Rl2g() 
    &   Kds_[fluid2_.name()+".structure"]
    &   structure_.Rg2l()
);

//- Correct turbulence
fluid1_.turbulence().correct();
fluid2_.turbulence().correct();

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
        -   fvm::Sp(fluid1_.contErr(), U1)
        +   fvm::Sp(dmdt12, U1)
        +   MRF_.DDt(alpha1*rho1, U1)
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
            
        +   fvm::Sp(Kd12, U1)
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
        -   fvm::Sp(fluid2_.contErr(), U2)
        +   fvm::Sp(dmdt21, U2)
        +   MRF_.DDt(alpha2*rho2, U2)
            //- The following is just a re-arrangement of div(nu*grad(U))
        +   (
                -   fvm::laplacian(alpha2*rho2*nuEff2, U2)
                -   fvc::div
                    (
                        alpha2*rho2*nuEff2&dev2(T(fvc::grad(U2))), 
                        "div(alphaRhoPhiNu,U)"
                    )
            )
        +   fvm::Sp(Kd12, U2)
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

    UEqn1.relax();
    fvOptions_.constrain(UEqn1);
    fvOptions_.correct(U1);
    UEqn2.relax();
    fvOptions_.constrain(UEqn2);
    fvOptions_.correct(U2);
}
