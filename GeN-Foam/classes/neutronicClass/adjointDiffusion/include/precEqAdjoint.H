forAll(prec_,precI)
{
    if(liquidFuel)
    {

    	prec_[precI] /= porosity_;//solve precs per volume of salt
        fvScalarMatrix precEq
        (
            fvm::ddt(porosity_*(1-eigenvalueNeutronics_),precStar_[precI])
            + fvm::Sp(lambda[precI]*porosity_, precStar_[precI])
            - adjointFlux_chiD_*lambda[precI]
            + fvm::div(-phi_, precStar_[precI])
            - fvm::laplacian(diffCoeffPrec_,precStar_[precI])                
        );
    	
        precEq.solve();
        prec_[precI] = precStar_[precI] * porosity_;//precs per total volume
    	prec_[precI] *= porosity_;//go back to precs per total volume
    }
    else
    {
        fvScalarMatrix precEq
        (
            fvm::ddt((1-eigenvalueNeutronics_),prec_[precI])
            + fvm::Sp(lambda[precI], prec_[precI])
            - adjointFlux_chiD_*lambda[precI]
        );
        precEq.solve();
    }
}

