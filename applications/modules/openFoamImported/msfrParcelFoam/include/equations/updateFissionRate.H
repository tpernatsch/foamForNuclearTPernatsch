{
    // Reset fields
    Rf_ *= 0.0;
    NuRf_ *= 0.0;
    EfRf_ *= 0.0;

    // Evaluate new fission rates based on current fluxes and cross-sections
    forAll(flux_, i)
    {
        Rf_ += Sf_[i]*flux_[i];
        NuRf_ += Nu_[i]*Sf_[i]*flux_[i];
        EfRf_ += Ef_[i]*Sf_[i]*flux_[i];
    }
}
