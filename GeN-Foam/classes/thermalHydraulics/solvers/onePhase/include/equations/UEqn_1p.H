//- The velocity predictor is missing in this implementation for
//  consistency with the two-phase treatment. It could be implemented
//  to improve convergence of specifc cases

//- Field refs
volVectorField& U(fluid_.U());
const volScalarField& alpha(fluid_);
volScalarField& rho(fluid_.rho());
fvVectorMatrix UEqn(U, dimMass*dimLength/dimTime/dimTime);
const volTensorField& Kds(FSPair_.Kd());

fluid_.turbulence().correct();

{
    volTensorField nuEff
    (
        structure_.tortuosity()
       *fluid_.turbulence().nuEff()()
    );

    UEqn =
    (
            fvm::div(fluid_.alphaRhoPhi(), U, "div(alphaRhoPhi,U)")
        +   fvm::SuSp(-fluid_.contErr(), U)
            //- The following is just a re-arrangement of div(nu*grad(U))
        +   (
                -   fvm::laplacian(alpha*rho*nuEff, U)
                -   fvc::div
                    (
                        alpha*rho*nuEff&dev2(T(fvc::grad(U))), 
                        "div(alphaRhoPhiNu,U)"
                    )
            )
            //- Implicit diagonal contribution of fluid-structure drag, 
            //  explicit off-diagonal contribution
        +   fvm::Sp((1.0/3.0)*tr(Kds), U) + (dev(Kds)&U)
        ==
            fvOptions_(alpha, rho, U)
    );

    if (momentumMode_ != momentumMode::faceCentered)
    {
        UEqn += fvm::ddt(alpha, rho, U);
    }

    UEqn.relax();
    fvOptions_.constrain(UEqn);
    fvOptions_.correct(U);
}
