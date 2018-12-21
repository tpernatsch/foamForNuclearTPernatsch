// This file sets constant values, that is, values that are not affected
// by feedbacks or at least are assumed to not change with time during
// the simulation.
 
forAll(flux,energyI)
{
    forAll(entries, zoneI)
    {
        label zone = zoneI;

        const word& name = entries[zoneI].keyword();

        Info << "Mesh zone: " << name  << endl;

        label zoneId = mesh_.cellZones().findZoneID(name);

        if(zoneId == -1)
        {
            Info << "Warning: mesh zone: " << name << " does not exists" << endl;
            Info << "Ignore this warning if it is a completely withdrawn control rod" << endl;
        }
        else
        {
            forAll(mesh_.cellZones()[zoneId], cellIlocal)
            {
                label cellIglobal = mesh_.cellZones()[zoneId][cellIlocal];
                zone = zoneI;

                //constant
                IV[energyI][cellIglobal] = IVList[zone][energyI];
                chiPrompt[energyI][cellIglobal] = chiPromptList[zone][energyI];
                chiDelayed[energyI][cellIglobal] = chiDelayedList[zone][energyI];
                BetaTot[cellIglobal] = BetaTotList[zone];

                fuelFraction[cellIglobal] = fuelFractionList[zone];

                
            }
        }
    }
    IV[energyI].correctBoundaryConditions();
    chiPrompt[energyI].correctBoundaryConditions();
    chiDelayed[energyI].correctBoundaryConditions();
    BetaTot.correctBoundaryConditions();
    fuelFraction.correctBoundaryConditions();
}
