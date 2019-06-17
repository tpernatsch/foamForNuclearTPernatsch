

forAll(prec_,precI)
{
    if(liquidFuel)
    {

        fvScalarMatrix precEq
        (
            fvm::ddt(porosity_*(1-eigenvalueNeutronics_),precStar_[precI])
            + fvm::Sp(lambda[precI]*porosity_, precStar_[precI])
            - neutroSource_/keff_*Beta[precI]
            + fvm::div(phi_, precStar_[precI])
            - fvm::laplacian(diffCoeffPrec_,precStar_[precI])                
        );
        precEq.solve();
        prec_[precI] = precStar_[precI] * porosity_;//precs per total volume
        prec_[precI].correctBoundaryConditions();

/*
	//correct but works on discontinuous fields
        fvScalarMatrix precEq
        (
            fvm::ddt((1-eigenvalueNeutronics_),prec_[precI])
            + fvm::Sp(lambda[precI], prec_[precI])
            - neutroSource_/keff_*Beta[precI]
            + fvm::div(phi_/fvc::interpolate(porosity_), prec_[precI])
            - fvm::laplacian(diffCoeffPrec_,prec_[precI])                
        );
        precEq.solve();
*/

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

