forAll(prec,precI)
{

    forAll(entries, zoneI)
    {
        const word& name = entries[zoneI].keyword();

        label zoneId = neutroMesh.cellZones().findZoneID(name);

        if(zoneId != -1)
        {

            const UList<label>& addr = neutroMesh.cellZones()[zoneId];
    
            IndirectList<scalar>(Beta[precI].internalField(), addr)
            = scalar(BetaList[zoneI][precI]);

            IndirectList<scalar>(lambda[precI].internalField(), addr)
            = scalar(lambdaList[zoneI][precI]);
        }

    }

    Beta[precI].correctBoundaryConditions();
    lambda[precI].correctBoundaryConditions();

}


