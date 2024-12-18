Info<< "Solving neutronics" << endl;

// Correct XS fields
xs_.correct(/*TFuel_, TClad_, rhoCool_, TCool_, disp_, TStructMech_*/);

#include "setXSFields.H"

// Read solver parameters from fvSolution
const dictionary& neutronTransport =
    mesh_.solutionDict().subDict("neutronTransport");

const bool integralPredictor =
    neutronTransport.lookupOrDefault("integralPredictor", false);
const bool implicitPredictor =
    neutronTransport.lookupOrDefault("implicitPredictor", false);
const bool aitkenAcceleration =
    neutronTransport.lookupOrDefault("aitkenAcceleration", false);
scalar neutronIterationResidual =
    neutronTransport.lookupOrDefault<scalar>
    (
        "neutronIterationResidual",
        0.00001
    );
int maxNeutronIterations =
    neutronTransport.lookupOrDefault<int>("maxNeutronIterations", 10);


// Correct flux with discontinuity factors
forAll(flux_, energyI)
{
    flux_[energyI] = fluxStar_[energyI] / discFactor[energyI];
}

// Compute total power
scalar pTot(1.0);

if
(
    eigenvalueNeutronics_ &&
    (mesh_.time().timeIndex() - mesh_.time().startTimeIndex()) == 1
)
{
    pTot = power();
    #include "normFluxes.H"
    pTot = power();
}

// Update mapped fields
if (liquidFuel_)
{
    phiPtr_() = fvc::flux(UPtr_() * alphaPtr_());
    diffCoeffPrecPtr_() = (alphatPtr_() + muPtr_() / xs_.ScNo()) / rhoCool_;
}
#include "initializeNeutroSource.H"
#include "precEq.H"
#include "initializeDelayedNeutroSource.H"

if (!eigenvalueNeutronics_)
{
    if (integralPredictor)
    {
        #include "integralPredictor.H"
    }
    if (implicitPredictor)
    {
        #include "implicitPredictor.H"
    }
}

int aitkenIterCounter = 0;
int aitkenIterNo = 0;
bool aitkenYesNo = true;
int neutroIter = 0;

#include "createFluxMatrices.H"

scalar neutroResidual = 0.0;

do
{
    Info<< "Neutron iteration no:  " << neutroIter << nl << endl;

    #include "fluxEq.H"

    #include "precEq.H"

    #include "initializeDelayedNeutroSource.H"

    neutroIter++;

    Info<< "neutroResidual (neutroIter: " << neutroIter << "): " << neutroResidual
        << endl;
        //<< " after Aitken no: " << aitkenIterNo << endl;
}
while
(
    (neutroResidual > neutronIterationResidual) &&
    (neutroIter < maxNeutronIterations) &&
    !eigenvalueNeutronics_
    // && !externalSourceNeutronics
);

// Reset and compute power density field
powerDensity_ *= 0.0;
forAll(flux_, energyI)
{
    powerDensity_ += flux_[energyI] * sigmaPow[energyI];
}

pTot = power();

Info<< "power: " << pTot << nl
    << "Final neutroResidual at neutroIter no "
    << neutroIter << ": "  << neutroResidual
    << " - pow:" << pTot
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

if (externalSourceNeutronics_)
{
    // Output the integrated power in reactor state
    reactorState_.set("pTarget", pTot);
}


residual_ = max(residual_, neutroResidual);

// Reset and compute one group field
oneGroupFlux_ *= 0.0;
forAll(flux_, energyI)
{
    oneGroupFlux_ += flux_[energyI];
}
