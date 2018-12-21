neutroResidual = 0.0;

if(aitkenAcceleration)
{
    #include "aitkenAcceleration.H"
}   

forAll(flux,energyI)
{
    Dalbedo = D[energyI];
    fluxStarAlbedo *= 0.0;
    #include "calcScatteringSource.H"
    neutroSource -= nuSigmaEff[energyI]*flux[energyI]*(1-eigenvalueNeutronics_);

    fvMatrixExt<scalar> neutroEq = neutroEqList[energyI];

    neutroEq = neutroEq
      - neutroSource*(1.0-BetaTot)*chiPrompt[energyI]/keff_
      - delayedNeutroSource*chiDelayed[energyI]
      - scatteringSourceExtra ;


    scalar initRes = neutroEq.initialResidual();
    neutroResidual = max(neutroResidual,initRes);

    Info<< "neutro Residual = " << initRes << nl << endl;

    neutroEq.solve();
   
    flux[energyI]=fluxStar[energyI]/discFactor[energyI];
    flux[energyI].correctBoundaryConditions();
    neutroSource += nuSigmaEff[energyI]*flux[energyI]*(1-eigenvalueNeutronics_);

}




