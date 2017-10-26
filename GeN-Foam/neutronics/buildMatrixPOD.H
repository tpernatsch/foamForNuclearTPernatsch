
    
{
    Info << "buildmatrixROM energy: " << energyI <<  endl;

    // Resize ROM matrices
    
    forAll(flux,energyJ)
    {
        maxPODBaseSize = max(maxPODBaseSize,base_fluxStar[energyJ].size());
    }
    forAll(flux,energyJ)//costretto a definire matrici della stessa dimensione
    {
        matrixROM[energyJ].setSize(maxPODBaseSize);
    }

    int PODBaseSize(base_fluxStar[energyI].size());
    for(int baseNoJ = 0; baseNoJ < PODBaseSize ; baseNoJ++ )
    {
        volScalarField Aphib = neutroEq.AmulWithBC(base_fluxStar[energyI][baseNoJ])();

        for(int baseNoI = 0; baseNoI < PODBaseSize ; baseNoI++ )
        {

            matrixROM[energyI][baseNoI][baseNoJ] = Foam::POD::projection(Aphib,base_fluxStar[energyI][baseNoI]);
        }
     
    }
    for(int baseNoI = PODBaseSize; baseNoI < maxPODBaseSize ; baseNoI++ )
    {
        if(matrixROM[energyI][baseNoI][baseNoI] == 0.0)
        {
            matrixROM[energyI][baseNoI][baseNoI] = 1.0;
        }
    }

}


