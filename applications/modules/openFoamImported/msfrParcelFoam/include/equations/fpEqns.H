{
    forAll(fp_, i)
    {
        // Stokes-Einstein molecular diffusivity (via FPDp_) + turbulent eddy diffusivity
        Dfp_[i] =
            constant::physicoChemical::k*T_
           /(3*constant::mathematical::pi*muc_*FPDp_[i])
          + turbulence_->nut()/Sct_;

        volScalarField& fpi = fp_[i];
        volScalarField& Dfpi = Dfp_[i];
        volScalarField& fpdi = fpd_[i];

        fvScalarMatrix fpEqni
        (
            fvm::ddt(fpi)
          + fvm::div(phi_, fpi)
          - fvm::laplacian(Dfpi, fpi)
          + fvm::Sp(FPLambda_[i], fpi)
         ==
            FPYield_[i]*Rf_
          + fvOptions_(fpi)
        );

        fpEqni.relax();
        fvOptions_.constrain(fpEqni);
        fpEqni.solve();
        fvOptions_.correct(fpi);

        #include "calcDeposition.H"
    }
}
