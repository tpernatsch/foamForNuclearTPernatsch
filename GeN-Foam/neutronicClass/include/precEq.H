forAll(prec_,precI)
{
    if(liquidFuel)
    {
        prec_[precI] /= porosity_;//solve precs per volume of salt
        fvScalarMatrix precEq
        (
            fvm::ddt(porosity_*(1-eigenvalueNeutronics_),prec_[precI])
            + fvm::Sp(lambda[precI]*porosity_, prec_[precI])
            - neutroSource_/keff_*Beta[precI]
            + fvm::div(phi_, prec_[precI])
            - fvm::laplacian(diffCoeffPrec_,prec_[precI])                
        );
        precEq.solve();
        prec_[precI] *= porosity_;//go back to precs per total volume
    }
    else
    {
        fvScalarMatrix precEq
        (
            fvm::ddt((1-eigenvalueNeutronics_),prec_[precI])
            + fvm::Sp(lambda[precI], prec_[precI])
            - neutroSource_/keff_*Beta[precI]
        );
        precEq.solve();
    }
}
