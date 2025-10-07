{
    fvScalarMatrix rhoEqn
    (
        fvm::ddt(rho_)
      + fvc::div(phi_)
      ==
        fvOptions_(rho_)
    );
 
    fvOptions_.constrain(rhoEqn);
 
    rhoEqn.solve();
 
    fvOptions_.correct(rho_);
}