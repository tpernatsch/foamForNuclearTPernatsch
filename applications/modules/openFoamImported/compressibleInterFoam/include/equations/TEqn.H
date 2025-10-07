{
    fvScalarMatrix TEqn
    (
        fvm::ddt(rho_, T_) + fvm::div(rhoPhi_, T_) - fvm::Sp(contErr, T_)
      - fvm::laplacian(turbulence_->alphaEff(), T_)
      + (
             divUp()// - contErr/rho_*p_
          +  (fvc::ddt(rho_, K_) + fvc::div(rhoPhi_, K_))() - contErr*K_
        )
       *(
           alpha1_()/mixture_.thermo1().Cv()()
         + alpha2_()/mixture_.thermo2().Cv()()
        )
     ==
        fvOptions_(rho_, T_)
    );

    TEqn.relax();

    fvOptions_.constrain(TEqn);

    TEqn.solve();

    fvOptions_.correct(T_);

    mixture_.correctThermo();
    mixture_.correct();
}
