//- Update reactivity coefficients


//- Check if you need to initialize precurosr distribution
bool initialize = false;
if(mesh_.time().timeIndex() == mesh_.time().startTimeIndex()+1)
{
    initialize = true;
}


#include "computeFeedbackFieldValues.H"
if (initialize)
{

        //- TFuelRef limited as it appears in a fraction denominator if doing
    //  fastNeutrons (for the Doppler coeff)

    TFuelRef_ =
       reactorState_.lookupOrDefault<scalar>("TFuelRef", TFuelValue);

    TCladRef_ =
        reactorState_.lookupOrDefault<scalar>("TCladRef", TCladValue);

    TCoolRef_ =
        reactorState_.lookupOrDefault<scalar>("TCoolRef", TCoolValue);

    rhoCoolRef_ =
        reactorState_.lookupOrDefault<scalar>("rhoCoolRef", rhoCoolValue);

    TStructRef_ =
        reactorState_.lookupOrDefault<scalar>("TStructRef", TStructValue);

    TStructMechRef_ =
        reactorState_.lookupOrDefault<scalar>("TStructMechRef", TStructMechValue);

    TDrivelineRef_ =
        reactorState_.lookupOrDefault<scalar>("TDrivelineRef", TDrivelineValue);

    //- Since a T*Ref_ are scalars, the value inside the dictionary is not
    //  updated at runTime automatically. The line below resets the SCALAR IN
    //  THE DICTIONARY at runTime so that the dictionary is written with the
    //  updated value.
    reactorState_.set("TFuelRef", TFuelRef_);
    reactorState_.set("TCladRef", TCladRef_);
    reactorState_.set("TCoolRef", TCoolRef_);
    reactorState_.set("rhoCoolRef", rhoCoolRef_);
    reactorState_.set("TStructRef", TStructRef_);
    reactorState_.set("TStructMechRef", TStructMechRef_);
    reactorState_.set("TDrivelineRef", TDrivelineRef_);
}

#include "correctReactivity.H"

//- Get current simulation time
const scalar timeOutput(mesh_.time().timeOutputValue());

//- Update external source power
if (externalSourceNeutronics_)
{
    if (externalSourceMode_ == "transient")
    {
        if (externalSourceModulationTimeProfile_.valid())
        {
            externalSourceModulation_ =
                externalSourceModulationTimeProfile_.value(timeOutput);
        }
    }
    else if (externalSourceMode_ == "powerMonitoring")
    {
        powerRecordPtr_.append(Foam::Tuple2<scalar, scalar>(t, power_));

        scalar powerDelayed(Foam::interpolationTable<scalar>::interpolateValue(
            powerRecordPtr_,
            timeOutput - powerModulationDelay_
        ));

        externalSourceModulation_ = externalSourceModulation_ * (
            1. + modulationFactor_ * (powerTarget_ - powerDelayed)/powerTarget_ * mesh_.time().deltaT().value()
        );
    }

    //- Update source strength
    externalSourcePower_ = externalSourcePowerRef_ * externalSourceModulation_;

    setBeamParameters();
}
else
{
    externalSourcePower_ = 0.0;
}

//- Update decayPower_, if available
// if (decayPowerPtr_.valid())
if (decayPowerTimeProfile_.valid())
{
    // decayPower_ = decayPowerPtr_->value
    // (
    //     mesh_.time().timeOutputValue() - decayPowerStartTime_
    // );
    decayPower_ = decayPowerTimeProfile_.value(timeOutput);
}

//- Update fissionPower_ and precursorPowers_ via actual pointKinetics model

//- Update old time step values if time step has changed
if (timeIndex_ != mesh_.time().timeIndex())
{
    timeIndex_ = mesh_.time().timeIndex();
    powerOld_ = power_;
    fissionPowerOld_ = fissionPower_;
    for (int i=0; i<delayedGroups_; i++)
    {
        precursorPowersOld_[i] = precursorPowers_[i];
    }
    externalSourcePowerOld_ = externalSourcePower_;
}
label n(delayedGroups_+1+1); // +1 for fissionPower; +1 for external source
scalar dt(mesh_.time().deltaT().value());

//- Coefficient matrix, source, solution field
SquareMatrix<scalar> A(n, 0.0);
List<scalar> B(n, 0.0);
List<scalar> x(n, 0.0);

//- First line of matrix has the coefficients of the discretized power
//  differential equation, the first element of the solution vector is the
//  power, the last element is the external source power
A[0][0] = 1.0/dt-(totalReactivity_-subcriticalIndex_-beta_)/promptGenerationTime_;
for (int i = 1; i < n-1; i++)
{
    A[0][i] = -lambdas_[i-1];
    B[0] = fissionPowerOld_/dt;
}
//- External Source
A[0][n-1] = -1/promptGenerationTime_;
A[n-1][n-1] = 1.0/dt;
B[n-1] = externalSourcePowerOld_/dt;

//- Rest of the matrix for the precursor coefficients
for (int i = 1; i < n-1; i++)
{
    A[i][0] = -betas_[i-1]/promptGenerationTime_;
    A[i][i] = 1.0/dt + lambdas_[i-1];
    B[i] = precursorPowersOld_[i-1]/dt;
}

//- Solve
solve(x, A, B);
fissionPower_ = x[0];
for (int i=1; i<n-1; i++)
{
    precursorPowers_[i-1] = x[i];
}

//- Update total power
power_ = fissionPower_ + decayPower_;

//- Re-scale fluxes, precursors by change in fissionPower
scalar fluxScaleFactor = fissionPower_/fissionPowerOld_;

for (int i = 0; i < energyGroups_; i++)
{
    fluxes_[i] *= fluxScaleFactor;
    fluxes_[i].correctBoundaryConditions();
}
oneGroupFlux_ *= fluxScaleFactor;
oneGroupFlux_.correctBoundaryConditions();

//- Re-scale powerDensity by change in total power
scalar powerDensityScaleFactor = power_/powerOld_;
powerDensity_ *= powerDensityScaleFactor;
powerDensity_.correctBoundaryConditions();
secondaryPowerDensity_ *= powerDensityScaleFactor;
secondaryPowerDensity_.correctBoundaryConditions();

//- Only possible to re-scale precursors if they were initially present,
//  so I'm using a forAllIter rather than a for, as delayedGroups_ is not
//  necessarily representative of the length of the prec_ PtrList
int i = 0;
forAllIter
(
    PtrList<volScalarField>,
    prec_,
    iter
)
{
    if (precursorPowersOld_[i] != 0)
    {
        volScalarField& precursori(iter());
        precursori *= precursorPowers_[i]/precursorPowersOld_[i];
    }
    i++;
}

//- Write precursorPowers to reactorState
if (TFuel_.mesh().time().write())
{
    reactorState_.set("power", power_);
    reactorState_.set("precursorPowers", precursorPowers_);
}

if(nuclearData_.found("updatePowerFields"))
{
    dictionary& subDict = nuclearData_.subDict("updatePowerFields");
    wordList regions = subDict.toc();
    forAll(regions, regI)
    {

        word fieldName = subDict.get<word>(regions[regI]);
        const dynamicFvMesh& regionMesh = mesh_.time().lookupObject<dynamicFvMesh>(regions[regI]);
        volScalarField& powerField = const_cast<volScalarField&>(regionMesh.lookupObject<volScalarField>(fieldName));
        powerField *= powerDensityScaleFactor;
    }
}

Info << endl << "pointKinetics: " << endl;
#include "pointKineticsInfo.H"
