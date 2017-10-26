pTot*=0;
forAll(flux,energyI)
{
    
    pTot += fvc::domainIntegrate(flux[energyI] * sigmaPow[energyI]).value();

}
Info << "pTot neutro: " << pTot << endl;

