

forAll(fluxStar,energyI)
{
    bool adjust = true;
    forAll(groupsWoDF, groupI)
    {
        if(energyI == groupsWoDF[groupI])
        {
            adjust = false;
            Info << "Not calculating DF for group: " << groupsWoDF[groupI] << endl;
        }
    }
    if(adjust)
    {
        volScalarField integrateYesNo
        (
            IOobject
            (
                "integrateYesNo",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("", dimensionSet(0,0,0,0,0,0,0),0.0),
            zeroGradientFvPatchScalarField::typeName
        );
        scalar totFluxSerpent = 0;
        scalar totFluxGeNFoam = 0;

        forAll(entries, zoneI)
        {
            
            const word& name = entries[zoneI].keyword();

            label zoneId = mesh_.cellZones().findZoneID(name);

            if(zoneId == -1)
            {

            }
            else
            {
                if(dfAdjustList[zoneI]==true)
                {
        
                    forAll(mesh_.cellZones()[zoneId], cellIlocal)
                    {
        
                        label cellIglobal = mesh_.cellZones()[zoneId][cellIlocal];
            
                        integrateYesNo[cellIglobal] = 1.0;
            
                    }
                    totFluxSerpent += integralFluxList[zoneI][energyI];
                
                }   
    
            }
        }

        totFluxGeNFoam = fvc::domainIntegrate(integrateYesNo*fluxStar[energyI]).value();

        integrateYesNo *= 0.0 ;
        forAll(entries, zoneI)
        {

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
                if(dfAdjustList[zoneI]==true)
                {
    
                    forAll(mesh_.cellZones()[zoneId], cellIlocal)
                    {
        
                        label cellIglobal = mesh_.cellZones()[zoneId][cellIlocal];
    
                        integrateYesNo[cellIglobal] = 1.0;
    
                    }

                    discFactorList[zoneI][energyI] =  (fvc::domainIntegrate(integrateYesNo*fluxStar[energyI]).value()/totFluxGeNFoam ) / (integralFluxList[zoneI][energyI]/totFluxSerpent);
                }
                
                Info << "discFactorList: zone " <<  entries[zoneI].keyword() << ", energy " << energyI << ": " << discFactorList[zoneI][energyI]  << endl;
    
                integrateYesNo *= 0.0 ;     
    
            }
        }
    }

}


