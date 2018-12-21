
scatteringSourceExtra *= 0.0;

forAll(flux,energyJ)
{
    if(energyJ != energyI)
    {
        scatteringSourceExtra += sigmaFromTo[energyJ][energyI] * flux[energyJ];
    }
}
