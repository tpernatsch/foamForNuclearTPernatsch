
PtrList<fvMatrixExt<scalar> > neutroEqList(energyGroups);

forAll(neutroEqList, energyI)
{
    Dalbedo = D[energyI];
    fluxStarAlbedo = fluxStar2[energyI]*(SP3Neutronics);
    neutroEqList.set
    (
        energyI,
        new fvMatrixExt<scalar>
        (
             fvm::ddt((1-eigenvalueNeutronics)*IV[energyI]/discFactor[energyI],fluxStar[energyI])
             - fvm::laplacian(D[energyI]/discFactor[energyI],fluxStar[energyI])
             - fvm::Sp(nuSigmaEff[energyI]/keff*(1.0-BetaTot)*chiPrompt[energyI]*(1-eigenvalueNeutronics)-sigmaDisapp[energyI]/discFactor[energyI],fluxStar[energyI])
        )
    );
}


