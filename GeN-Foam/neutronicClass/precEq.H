forAll(prec,precI)
{

    prec[precI] /= porosity;//solve precs per volume of salt
        fvScalarMatrix precEq
        (
            fvm::ddt(porosity*(1-eigenvalueNeutronics_),prec[precI])
            + fvm::Sp(lambda[precI]*porosity, prec[precI])
            - neutroSource/keff_*Beta[precI]
           // + fvm::div(phi, prec[precI])
           // - fvm::laplacian(diffCoeffPrec,prec[precI])                
        );

    precEq.solve();
    prec[precI] *= porosity;//go back to precs per total volume

}
