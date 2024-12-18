Info<< "Solving neutronics" << endl;

// Correct XS fields
xs_.correct(/*TFuel_, TClad_, rhoCool_, TCool_, disp_, TStructMech_*/);

#include "setXSFields.H"

// Read solver parameters from fvSolution
const dictionary& neutronTransport = mesh_.solutionDict().subDict("neutronTransport");

const bool integralPredictor = neutronTransport.lookupOrDefault("integralPredictor", false);
const bool implicitPredictor = neutronTransport.lookupOrDefault("implicitPredictor", false);
const bool aitkenAcceleration = neutronTransport.lookupOrDefault("aitkenAcceleration", false);
scalar neutronIterationResidual =
    neutronTransport.lookupOrDefault<scalar>("neutronIterationResidual", 0.00001);
int maxNeutronIterations =
    neutronTransport.lookupOrDefault<int>("maxNeutronIterations", 10);

// For Aitken
List<scalar> integralFlux0(xs_.energyGroups());
List<scalar> integralFlux1(xs_.energyGroups());
List<scalar> integralFlux2(xs_.energyGroups());

forAll(flux_, energyI)
{
    flux_[energyI] = (fluxStar_[energyI] - fluxStar2_[energyI] * 2.0) / discFactor[energyI];
}

scalar pTot(1.0);

if
(
    eigenvalueNeutronics_ &&
    (mesh_.time().timeIndex() - mesh_.time().startTimeIndex()) == 1
)
{
    pTot = power();
    #include "normFluxesSP3.H"
    pTot = power();
}

//- Update mapped fields
if (liquidFuel_)
{
    phiPtr_() = fvc::flux(UPtr_() * alphaPtr_());
    diffCoeffPrecPtr_() = (alphaPtr_() + muPtr_() / xs_.ScNo()) / rhoCool_;
}

#include "initializeNeutroSource.H"
#include "precEq.H"

#include "initializeDelayedNeutroSource.H"

if (!eigenvalueNeutronics_)
{
    if (integralPredictor)
    {
        #include "integralPredictorSP3.H"
    }
    if (implicitPredictor)
    {
        #include "implicitPredictorSP3.H"
    }
}

int aitkenIterCounter = 0;
int aitkenIterNo = 0;
bool aitkenYesNo = true;

int neutroIter = 0;

#include "createFluxMatricesSP3.H"

scalar neutroResidual = 0.0;

do
{
    Info << "Neutron iteration no:  " << neutroIter << nl << endl;

    #include "fluxEqSP3.H"
    #include "precEq.H"
    #include "initializeDelayedNeutroSource.H"

    neutroIter++;

    Info<< "neutroResidual (neutroIter: " << neutroIter << "): "
        << neutroResidual << " after Aitken no: " << aitkenIterNo
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

Info<< "Final neutroResidual at neutroIter no " << neutroIter << ": "
    << neutroResidual << " - pow:" << pTot
    << endl;

//- Determine power densities
secondaryPowerDensity_ = powerDensity_  / max(secondaryPowerVolumeFraction, SMALL) * fractionToSecondaryPower  ;
powerDensity_ /= max(fuelFraction, SMALL) * (1.0 - fractionToSecondaryPower);


if (eigenvalueNeutronics_)
{
    #include "calcKeff.H"

    #include "normFluxesSP3.H"
    xs_.adjustDiscFactors(fluxStar_);
}

residual_ = max(residual_,neutroResidual);

oneGroupFlux_ *= 0.0;
forAll(flux_, energyI)
{
    oneGroupFlux_+=flux_[energyI];
}
