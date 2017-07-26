{
    thermo.correct();
    volScalarField& he = thermo.he();

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho*porosity, he) + fvm::div(phi, he)
      + fvc::ddt(rho*porosity, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rho), U),
                p,
                "div(phiv,p)"
            )
          : -dpdt
        )
      - fvm::laplacian(turb.alphaEff(), he)
     ==
          porousMedium.externalVolHeatSource()
        - pow((pow(porousMedium.externalHeatTransferCoefficient(),-1)+pow(porousMedium.heatTransferCoefficient(),-1)),-1)* porousMedium.volumetricAreaSS()*(thermo.T()-porousMedium.externalT())
     	+ subscaleFuel.heatSources(thermo.T(),porousMedium.heatTransferCoefficient(),volFuelPower) * porousMedium.volumetricAreaFuel()*(1-liquidFuel)
	+ liquidFuel*volFuelPower*porosity
     	- porousMedium.heatTransferCoefficient()* porousMedium.volumetricAreaSS()*(thermo.T()-Tstructures)

    );

    EEqn.relax();

    scalar EEqnResidual(EEqn.solve(fluidMesh.solver(he.select(finalIter))).initialResidual());
    residual = max(residual,EEqnResidual);

    thermo.correct();

    if(!compr)
    {
    	rhok = 1.0 - beta*(thermo.T() - Tref);
    }

    Info<< "Min/max T:" << min(thermo.T()).value() << ' '
    << max(thermo.T()).value() << endl;

     TinnerFuel = subscaleFuel.fuelIn();
     TouterFuel = subscaleFuel.fuelOut();
     TinnerClad = subscaleFuel.cladIn();
     TouterClad = subscaleFuel.cladOut();

     //make avarege fuel temps outside of the fuel equals to the average in the fuel
     volScalarField fuelYesNo = TinnerFuel/max(TinnerFuel,dimensionedScalar("", dimensionSet(0,0,0,1,0,0,0),SMALL));
     TavFuel = fvc::domainIntegrate(TavFuel*fuelYesNo)/max(fvc::domainIntegrate(TavFuel*fuelYesNo/max(TavFuel*fuelYesNo,dimensionedScalar("", dimensionSet(0,0,0,1,0,0,0),SMALL))),dimensionedScalar("", dimensionSet(0,3,0,0,0,0,0),SMALL));
     TavClad = fvc::domainIntegrate(TavClad*fuelYesNo)/max(fvc::domainIntegrate(TavClad*fuelYesNo/max(TavClad*fuelYesNo,dimensionedScalar("", dimensionSet(0,0,0,1,0,0,0),SMALL))),dimensionedScalar("", dimensionSet(0,3,0,0,0,0,0),SMALL));
     

     subscaleFuel.fuelAverage(TavFuel);
     subscaleFuel.cladAverage(TavClad);


    //solve for Temp of structures
    fvScalarMatrix structureEqn
    (
        fvm::ddt(porousMedium.externalRhoCp(), Tstructures)
        ==
        porousMedium.heatTransferCoefficient()* porousMedium.volumetricAreaSS()*(thermo.T()-Tstructures)
    );

    structureEqn.solve();

}


