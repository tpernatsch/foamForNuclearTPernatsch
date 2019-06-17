forAll(prec_,precI)
{
    if(liquidFuel)
    {
/*
        prec_[precI] /= porosity_;//solve precs per volume of salt
        prec_[precI].correctBoundaryConditions();
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
        prec_[precI].correctBoundaryConditions();
*/

	//Problem with above is that ddt will use the old value, which is
	//not divided by porosity
	//quick fix here below: correct but works on discontinuous fields
        fvScalarMatrix precEq
        (
            fvm::ddt((1-eigenvalueNeutronics_),prec_[precI])
            + fvm::Sp(lambda[precI], prec_[precI])
            - neutroSource_/keff_*Beta[precI]
            + fvm::div(phi_/fvc::interpolate(porosity_), prec_[precI])
            - fvm::laplacian(diffCoeffPrec_,prec_[precI])                
        );
        precEq.solve();


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

