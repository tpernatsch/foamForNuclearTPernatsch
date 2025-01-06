    fvVectorMatrix UEqn
    (
        fvm::ddt(rho_, U_) + fvm::div(rhoPhi_, U_)
      - fvm::Sp(contErr, U_)
      + MRF_.DDt(rho_, U_)
      + turbulence_->divDevRhoReff(U_)
     ==
        fvOptions_(rho_, U_)
    );

    UEqn.relax();

    fvOptions_.constrain(UEqn);

    if (pimple_.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    mixture_.surfaceTensionForce()
                  - ghf_*fvc::snGrad(rho_)
                  - fvc::snGrad(p_rgh_)
                ) * mesh_.magSf()
            )
        );

        fvOptions_.correct(U_);

        K_ = 0.5*magSqr(U_);
    }
