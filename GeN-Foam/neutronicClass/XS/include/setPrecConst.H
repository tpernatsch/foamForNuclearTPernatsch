forAll(Beta_,precI)
{

    forAll(entries, zoneI)
    {
        const word& name = entries[zoneI].keyword();

        label zoneId = mesh_.cellZones().findZoneID(name);

        if(zoneId != -1)
        {

            const UList<label>& addr = mesh_.cellZones()[zoneId];
    
            IndirectList<scalar>(Beta_[precI].internalField(), addr)
            = scalar(BetaList[zoneI][precI]);

            IndirectList<scalar>(lambda_[precI].internalField(), addr)
            = scalar(lambdaList[zoneI][precI]);
        }

    }

    Beta_[precI].correctBoundaryConditions();
    lambda_[precI].correctBoundaryConditions();

}


