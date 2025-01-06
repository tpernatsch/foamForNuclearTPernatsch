Info<< "Solving neutronics" << endl;

// Correct XS fields
xs_.correct(TFuel_, TClad_, rhoCool_, TCool_, disp_, TStructMech_);

#include "setXSFields.H"

// Read solver parameters from fvSolution
const dictionary& neutronTransport = mesh_.solutionDict().subDict("neutronTransport");

scalar neutronIterationResidual =
    neutronTransport.lookupOrDefault<scalar>("neutronIterationResidual", 0.00001);
int maxNeutronIterations =
    neutronTransport.lookupOrDefault<int>("maxNeutronIterations", 10);

scalar pTot(1.0);

forAll(flux_, energyI)
{
    flux_[energyI] = fluxStar_[energyI] / discFactor[energyI];
}

if
(
    eigenvalueNeutronics_ &&
    (mesh_.time().timeIndex() - mesh_.time().startTimeIndex()) == 1
)
{
    pTot = power();
    #include "normFluxesAdjoint.H"
    pTot = power();
}

#include "initializeNeutroSourceAdjoint.H"
#include "initializeDelayedNeutroSourceAdjoint.H"

#include "initialize_adjointFlux_chiD.H"

int neutroIter = 0;

#include "createFluxMatricesAdjoint.H"

scalar neutroResidual = 0.0;

//-Update mapped fields
if (liquidFuel_)
{
    phiPtr_() = fvc::flux(UPtr_() * alphaPtr_());
    diffCoeffPrecPtr_() = (alphaPtr_() + muPtr_() / xs_.ScNo()) / rhoCool_;
}


do
{
    Info << "Neutron iteration no:  " << neutroIter << nl << endl;

    #include "fluxEqAdjoint.H"
    #include "precEqAdjoint.H"
    #include "initialize_adjointFlux_chiD.H"
    #include "initializeDelayedNeutroSourceAdjoint.H"

    neutroIter++;

    Info<< "neutroResidual (neutroIter: " << neutroIter << "): " << neutroResidual
        << endl;
}
while
(
    (neutroResidual > neutronIterationResidual) &&
    (neutroIter < maxNeutronIterations) &&
    !eigenvalueNeutronics_
);

powerDensity_ *= 0.0;
forAll(flux_, energyI)
{
    powerDensity_ += flux_[energyI] * sigmaPow[energyI];
}

pTot = power();

Info<< "power: " << pTot << endl;

Info<< "Final neutroResidual at neutroIter  no " << neutroIter << ": "
    << neutroResidual << " - pow:" << pTot
    << endl;

// Determine power densities
secondaryPowerDensity_ = powerDensity_ / max(secondaryPowerVolumeFraction, SMALL) * fractionToSecondaryPower;
powerDensity_ /= max(fuelFraction, SMALL) * (1.0 - fractionToSecondaryPower);


if (eigenvalueNeutronics_)
{
    #include "calcKeff.H"

    #include "normFluxes.H"
    xs_.adjustDiscFactors(fluxStar_);
}

residual_ = max(residual_,neutroResidual);

oneGroupFlux_ *= 0.0;
forAll(flux_,energyI)
{
    oneGroupFlux_ += flux_[energyI];
}

#include "initializeNeutroSourceAdjoint.H"
#include "initialize_adjointFlux_chiD.H"
#include "initializeDelayedNeutroSourceAdjoint.H"
