neutroResidual = 0.0;
/*
if(aitkenAcceleration)
{
    #include "aitkenAcceleration.H"
}   
*/



forAll(flux_,energyI)
{
    #include "calcSourceSN.H"
    scalar innerResidual = 0.0;
    int innerIter = 0;

    do
    {
        innerResidual = 0.0;
        forAll(angularFlux_[energyI], directionI)
        {
            #include "calcThisIntraSource.H"
            facePhi_ = facePhis_[directionI];
            scalar localResidual =  solve
            (
                //neutroEqList[energyI][directionI]
                fvm::ddt((1-eigenvalueNeutronics_)*IV[energyI],angularFlux_[energyI][directionI])
                + fvm::div(facePhi_, angularFlux_[energyI][directionI],"div(facePhi_,angularFlux_)")
                + fvm::Sp(sigmaDisapp[energyI] + sigmaFromToP[0][energyI][energyI], angularFlux_[energyI][directionI])
                ==
                extraSources[directionI]
                + intraSources[directionI]
                + delayedNeutroSource_*chiDelayed[energyI]
            ).initialResidual();
            innerResidual = max(innerResidual, localResidual);
            neutroResidual = max(neutroResidual, localResidual);
            if(eigenvalueNeutronics_) angularFlux_[energyI][directionI] = max(angularFlux_[energyI][directionI],(angularFlux_[energyI][directionI]*0.0));
            
        }
        #include "updateThisScalarFlux.H"
        innerIter++;
    }
    while((innerResidual > (neutroResidual*innToOutTolerance))&&(innerIter < maxInnerIter));
    intraSources.clear();
    extraSources.clear();

}
/*
forAll(flux_,energyJ)
{
    if(energyJ != energyI)
    {
        scatteringSourceExtra_ += sigmaFromTo[energyJ][energyI] * flux_[energyJ];
    }
}
neutroSource_ *= 0.0;
forAll(flux_,energyI)
{
    neutroSource_ += nuSigmaEff[energyI]*flux_[energyI];
}
forAll(neutroEqList, energyI)
{

    Dalbedo_ = D[energyI];
    fluxStarAlbedo_ *= 0.0;
    neutroEqList.set
    (
        energyI,
        new fvMatrixExt<scalar>
        (
             fvm::ddt((1-eigenvalueNeutronics_)*IV[energyI]/discFactor[energyI],fluxStar_[energyI])
             - fvm::laplacian(D[energyI]/discFactor[energyI],fluxStar_[energyI])
             - fvm::Sp(nuSigmaEff[energyI]/keff_*(1.0-BetaTot)*chiPrompt[energyI]*(1-eigenvalueNeutronics_)-sigmaDisapp[energyI]/discFactor[energyI],fluxStar_[energyI])
        )
    );
}
  neutroSource_ -= nuSigmaEff[energyI]*flux_[energyI]*(1-eigenvalueNeutronics_);

    fvMatrixExt<scalar> neutroEq  = neutroEqList[energyI]
      - neutroSource_*(1.0-BetaTot)*chiPrompt[energyI]/keff_
      - delayedNeutroSource_*chiDelayed[energyI]
      - scatteringSourceExtra_ ;
      */
