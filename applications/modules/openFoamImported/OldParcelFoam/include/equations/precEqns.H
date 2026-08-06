{
    volScalarField Dprec("Dprec", turbulence_->nu()/Sc_ + turbulence_->nut()/Sct_);

    forAll(prec_, i)
    {
        volScalarField& preci = prec_[i];

        fvScalarMatrix precEqni
        (
            fvm::ddt(preci)
          + fvm::div(phi_, preci)
          - fvm::laplacian(Dprec, preci)
          + fvm::Sp(precLambda_[i], preci)
         ==
            precBeta_[i]*NuRf_
          + fvOptions_(preci)
        );

        precEqni.relax();
        fvOptions_.constrain(precEqni);
        precEqni.solve();
    }
}
