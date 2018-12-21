
PtrList<fvMatrixExt<scalar> > neutroEqList(energyGroups);

forAll(neutroEqList, energyI)
{
    Dalbedo = D[energyI];
    fluxStarAlbedo *= 0.0;
    neutroEqList.set
    (
        energyI,
        new fvMatrixExt<scalar>
        (
             fvm::ddt((1-eigenvalueNeutronics_)*IV[energyI]/discFactor[energyI],fluxStar[energyI])
             - fvm::laplacian(D[energyI]/discFactor[energyI],fluxStar[energyI])
             - fvm::Sp(nuSigmaEff[energyI]/keff_*(1.0-BetaTot)*chiPrompt[energyI]*(1-eigenvalueNeutronics_)-sigmaDisapp[energyI]/discFactor[energyI],fluxStar[energyI])
        )
    );
}


