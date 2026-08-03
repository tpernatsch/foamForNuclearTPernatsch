{
    volScalarField Ddec("Ddec", turbulence_->nu()/Sc_ + turbulence_->nut()/Sct_);

    forAll(dec_, i)
    {
        volScalarField& deci = dec_[i];

        fvScalarMatrix decEqni
        (
            fvm::ddt(deci)
          + fvm::div(phi_, deci)
          - fvm::laplacian(Ddec, deci)
          + fvm::Sp(decLambda_[i], deci)
         ==
            decBeta_[i]*EfRf_
          + fvOptions_(deci)
        );

        decEqni.relax();
        fvOptions_.constrain(decEqni);
        decEqni.solve();
    }
}
