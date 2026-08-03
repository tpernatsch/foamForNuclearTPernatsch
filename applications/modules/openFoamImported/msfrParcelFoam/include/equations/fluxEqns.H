{
    forAll(flux_, i)
    {
        volScalarField& fluxi = flux_[i];

        *diffCoeffAlbedo_ = D_[i];

        fluxi.correctBoundaryConditions();

        fvScalarMatrix fluxEqni
        (
            fvm::ddt(invVel_[i], fluxi)
          - fvm::laplacian(D_[i], fluxi)
          - fvm::Sp(M_[i], fluxi)
         ==
            neutronSource_[i]
        );

        fluxEqni.relax();
        fvOptions_.constrain(fluxEqni);
        fluxEqni.solve();
    }

    // Update total flux
    fluxTot_ *= 0.0;
    forAll(flux_, i)
    {
        fluxTot_ += flux_[i];
    }
}
