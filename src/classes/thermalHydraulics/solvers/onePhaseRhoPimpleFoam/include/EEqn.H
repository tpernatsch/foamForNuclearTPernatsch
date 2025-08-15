{
    volScalarField& he = thermo_.he();

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho_, he) + fvm::div(phi_, he)
      + fvc::ddt(rho_, K_) + fvc::div(phi_, K_)
      + (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi_/fvc::interpolate(rho_), U_),
                p_,
                "div(phiv,p)"
            )
          : -dpdt_
        )
      - fvm::laplacian(turbulence_->alphaEff(), he)
     ==
        fvOptions_(rho_, he)
    );

    if (MRF_.active())
    {
        EEqn += fvc::div(MRF_.phi(), p_);
    }

    EEqn.relax();

    fvOptions_.constrain(EEqn);

    EEqn.solve();

    fvOptions_.correct(he);

    thermo_.correct();

    kappaEff_ = 
    (
      thermo_.kappaEff(turbulence_->alphat())
    );
}