
const volScalarField& porosity = porousMedium.voidFraction();
realU = U / porosity;
realU.correctBoundaryConditions();
volScalarField& p = thermo.p();

if(compr)
{
    rhok = 1.0;
}
else
{
    rhok = 1.0 - beta*(thermo.T() - Tref);
}

//make avarege fuel temps outside of the fuel equals to the average in the fuel

tmp<volScalarField> tFuelYesNo = TinnerFuel/max(TinnerFuel,dimensionedScalar("", dimensionSet(0,0,0,1,0,0,0),SMALL));
const volScalarField& fuelYesNo = tFuelYesNo();
TavFuel = fvc::domainIntegrate(TavFuel*fuelYesNo)/max(fvc::domainIntegrate(TavFuel*fuelYesNo/max(TavFuel*fuelYesNo,dimensionedScalar("", dimensionSet(0,0,0,1,0,0,0),SMALL))),dimensionedScalar("", dimensionSet(0,3,0,0,0,0,0),SMALL));
TavClad = fvc::domainIntegrate(TavClad*fuelYesNo)/max(fvc::domainIntegrate(TavClad*fuelYesNo/max(TavClad*fuelYesNo,dimensionedScalar("", dimensionSet(0,0,0,1,0,0,0),SMALL))),dimensionedScalar("", dimensionSet(0,3,0,0,0,0,0),SMALL));
subscaleFuel.fuelAverage(TavFuel);
subscaleFuel.cladAverage(TavClad);
tFuelYesNo.clear();


const dimensionedScalar initialMass
(
    "initialMass",
    dimMass,
    fvc::domainIntegrate(rho).value()
);

//For incompressible
label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    fluidMesh.solutionDict().subDict("PIMPLE"),
    pRefCell,
    pRefValue
);
if(!compr)
{
    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
    }
}

bool& resReachedFluid = residualReachedFluid;
bool& residualControlUsed = residualControlUsedFluid;
/*
bool& resReachedFluid = residualReachedFluid[i];
bool& residualControlUsed = residualControlUsedFluid[i];
*/
