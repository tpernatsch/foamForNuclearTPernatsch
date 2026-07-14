//- Get solution control
const dictionary& neutronTransport =
mesh_.solutionDict().subDict("neutronTransport");

scalar neutronIterationResidual =
neutronTransport.lookupOrDefault<scalar>
(
    "neutronIterationResidual",
    0.00001
);
int maxNeutronIterations =
neutronTransport.lookupOrDefault<int>("maxNeutronIterations", 10);

//- Update reactivity coefficients
#include "computeFeedbackFieldValues.H"



//- Update old time step values if time step has changed
if(timeIndex_ != mesh_.time().timeIndex())
{
    timeIndex_ = mesh_.time().timeIndex();
    powerOld_ = power_;
}

scalar dt(mesh_.time().deltaT().value());

//- prepare fields necessary to weight precursors
PtrList<scalar> weightedPrec(delayedGroups_);


//- Iterate precursors and power till sufficiently small prec residual
int neutroIter = 0;
scalar neutroResidual = 0.0;

//- Reference feedback state (reference temperatures/densities and the
//  precursor-equilibrium reactivity) is established/loaded at the first step
//  of every run, independently of initPrecursorsLiquidFuel_. This mirrors the
//  solid-fuel path (solvePointKinetics.H) and, together with the reads below,
//  allows a restart to recover these cheap scalars from reactorState without
//  repeating the expensive equilibrium precursor solve.
bool initReferenceState =
    (mesh_.time().timeIndex() == mesh_.time().startTimeIndex()+1);

//- Whether to (re)compute the equilibrium precursor distribution via the
//  pseudo time-stepping solve below. This is the expensive step (the ~hours
//  process) and is only performed on explicit request.
bool initialize =
    initPrecursorsLiquidFuel_ && initReferenceState && neutroIter==0;


if (initReferenceState)
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

    //- Load the precursor-equilibrium reactivity persisted by a previous
    //  initialization run. When it is present (a restart from a run that
    //  computed it) it is read back here, so the transient can start from the
    //  correct steady state without redoing the equilibrium precursor solve.
    //  When it is absent and we are not about to compute it now, the transient
    //  would start from a non-critical state: warn the user.
    if (reactorState_.found("precEquilibriumReactivity"))
    {
        precEquilibriumReactivity_ =
            reactorState_.get<scalar>("precEquilibriumReactivity");
    }
    else if (!initialize)
    {
        WarningInFunction
            << "precEquilibriumReactivity not found in reactorState while "
            << "initPrecursorsLiquidFuel is false. The precursor equilibrium "
            << "state has never been computed, so this liquid-fuel transient "
            << "will start from a non-equilibrium (non-critical) condition. "
            << "Run once with initPrecursorsLiquidFuel = true to generate it, "
            << "then restart with the flag off to reuse it." << endl;
    }
}

#include "correctReactivity.H"


//- Correct the diffusion coefficient and flux
diffCoeffPrecPtr_() = (alphatPtr_()+muPtr_()/ScNo_)/rhoCool_; // exactly the same
phiPtr_()=fvc::flux(UPtr_()*alphaPtr_());

//Calc equilibrium precPK at the steady-state using pseudo time stepping
if(initialize)
{
        //- Determine maximum number of iterations among precursors
        label maxPseudo(0);
        for (label i = 0; i < pseudoTimeSteps_.size(); ++i)
        {
            maxPseudo = max(maxPseudo, pseudoTimeSteps_[i]);
        }

        if (maxPseudo > 0)
        {
            scalarList finalPseudoResiduals(delayedGroups_, -1.0);

            Info << "Pseudo time stepping initialization for precursors:" << nl
                 << "  Max iterations: " << maxPseudo << nl
                 << "  Per-group iterations: ";
            forAll(pseudoTimeSteps_, i) Info << pseudoTimeSteps_[i] << " ";
            Info << nl << "  Per-group dt: ";
            forAll(pseudoTimeStepSizes_, i) Info << pseudoTimeStepSizes_[i] << " ";
            Info << nl << "  Per-group relaxation: ";
            forAll(pseudoTimeRelaxations_, i) Info << pseudoTimeRelaxations_[i] << " ";
            Info << endl;

            //- Store old values for relaxation
            PtrList<volScalarField> precPKStarOld(delayedGroups_);
            forAll(precPK_, precI)
            {
                precPKStarOld.set(precI, new volScalarField(precPKStar_[precI]));
            }

            //- Pseudo time stepping loop up to max across groups. Each group
            //  only participates for its configured number of iterations
            for (int pseudoIter = 0; pseudoIter < maxPseudo; ++pseudoIter)
            {
                forAll(precPK_, precI)
                {
                    if (pseudoIter >= pseudoTimeSteps_[precI])
                    {
                        continue;
                    }

                    scalar pseudoDt = pseudoTimeStepSizes_[precI];

                    fvMatrix<scalar> precEq
                    (
                        fvm::Sp(
                            dimensionedScalar("", dimless/dimTime, 1.0/pseudoDt)*alphaPtr_()
                            + dimensionedScalar("", dimless/dimTime, lambdas_[precI])*alphaPtr_(),
                            precPKStar_[precI]
                        )
                        - isFuel_*initOneGroupFluxN_*dimensionedScalar("", dimPower, power_)*betas_[precI]/
                          dimensionedScalar("", dimTime, promptGenerationTime_)
                        + fvm::div(phiPtr_(), precPKStar_[precI], "div(phi_,precStar_)")
                        - fvm::laplacian(diffCoeffPrecPtr_(), precPKStar_[precI])
                    );

                    //- Add the source term from previous iteration: (1/pseudoDt)*alpha*precPKStarOld
                    precEq.source() +=
                        (1.0/pseudoDt) * alphaPtr_().internalField() * precPKStarOld[precI].internalField()
                        * precPKStar_[precI].mesh().V();

                    precEq.solve();

                    //- Apply per-group relaxation if specified
                    if (pseudoTimeRelaxations_[precI] < 1.0)
                    {
                        precPKStar_[precI] =
                            pseudoTimeRelaxations_[precI] * precPKStar_[precI]
                            + (1.0 - pseudoTimeRelaxations_[precI]) * precPKStarOld[precI];
                    }

                    //- Store per-group ending residual to guide pseudo-time tuning
                    scalar pseudoDenom =
                        fvc::domainIntegrate(mag(precPKStarOld[precI])).value()
                        + SMALL;
                    scalar pseudoNumer =
                        fvc::domainIntegrate
                        (
                            mag(precPKStar_[precI] - precPKStarOld[precI])
                        ).value();
                    finalPseudoResiduals[precI] = pseudoNumer/pseudoDenom;

                    //- Update old for next pseudo iteration
                    precPKStarOld[precI] = precPKStar_[precI];
                }

                if (pseudoIter % max(1, maxPseudo/10) == 0 || pseudoIter == maxPseudo - 1)
                {
                    Info << "  Pseudo time step " << (pseudoIter + 1)
                         << " / " << maxPseudo << endl;
                }
            }

            Info << "Initialization complete." << nl
                 << "  Final residuals by precursor group:" << nl;
            forAll(finalPseudoResiduals, i)
            {
                if (finalPseudoResiduals[i] < 0.0)
                {
                    Info << "    group " << i << ": not solved (pseudoTimeSteps = 0)" << nl;
                }
                else
                {
                    Info << "    group " << i << ": " << finalPseudoResiduals[i] << nl;
                }
            }
            Info << endl;
        }
    else
    {
        //- If pseudoTimeSteps = 0, do steady state
        Info << "Using standard steady-state solver for precursor initialization" << nl
             << endl;

        forAll(precPK_, precI)
        {
            fvMatrix<scalar> precEq
            (
                fvm::Sp(dimensionedScalar("", dimless/dimTime, lambdas_[precI])*alphaPtr_(), precPKStar_[precI])
                - isFuel_*initOneGroupFluxN_*dimensionedScalar("", dimPower, power_)*betas_[precI]/
                  dimensionedScalar("", dimTime, promptGenerationTime_)
                + fvm::div(phiPtr_(), precPKStar_[precI], "div(phi_,precStar_)")
                - fvm::laplacian(diffCoeffPrecPtr_(), precPKStar_[precI])
            );
            precEq.solve();
            precPKStar_[precI].storePrevIter();
        }
    }
}

/*
//Calc the precEquilibriumReactivity
if(initialize)
{
    precEquilibriumReactivity_= 0.0;
    precEquilibriumReactivity_ += beta_;

    forAll(precPKStar_,precI)
    {
        //- Conver to precs per total volume
        precPK_[precI] = precPKStar_[precI] * alphaPtr_();
        precPK_[precI].correctBoundaryConditions();

        weightedPrec.set(precI,new scalar(fvc::domainIntegrate(
        initOneGroupFluxN_*precPK_[precI]).value()/domainIntegratedInitOneGroupFluxN_ ));

        precEquilibriumReactivity_ -= weightedPrec[precI] * lambdas_[precI] *
        promptGenerationTime_ / powerOld_;
    }
    #include "correctReactivity.H"
    Info << "precEquilibriumReactivity[pcm] =  " << precEquilibriumReactivity_ *1e+5 << nl << endl;
}
*/



//Calc beta for liquid fuel
liquidFuelBeta_ *= 0.0;
forAll(precPKStar_,precI)
{
    //- Convert to precs per total volume
    precPK_[precI] = precPKStar_[precI] * alphaPtr_();
    precPK_[precI].correctBoundaryConditions();

    weightedPrec.set(precI,new scalar(fvc::domainIntegrate(
    initOneGroupFluxN_*precPK_[precI]).value()/domainIntegratedInitOneGroupFluxN_ ));

    liquidFuelBeta_ += weightedPrec[precI] * lambdas_[precI] *
    promptGenerationTime_ / powerOld_;
}

//Calc the precEquilibriumReactivity
if(initialize)
{
    precEquilibriumReactivity_ = 0.0;
    precEquilibriumReactivity_ += beta_;
    precEquilibriumReactivity_ -= liquidFuelBeta_;

    //- Persist it in reactorState (AUTO_WRITE) so a later restart can reload
    //  it (see the reference-state block above) instead of recomputing the
    //  equilibrium precursor distribution.
    reactorState_.set("precEquilibriumReactivity", precEquilibriumReactivity_);

    #include "correctReactivity.H"
}

Info << "Beta[pcm] =  " << beta_ *1e+5 << nl << endl;
Info << "liquidFuelBeta[pcm] =  " << liquidFuelBeta_ *1e+5 << nl << endl;
Info << "precEquilibriumReactivity[pcm] =  " << precEquilibriumReactivity_ *1e+5 << nl << endl;


do
{
    Info << "Neutron iteration no:  " << neutroIter << nl << endl;

    //- Solve for precursors
    neutroResidual = 0.0;


    forAll(precPK_,precI)
    {
        //- Solve precursors
        fvMatrix<scalar> precEq
        (
            fvm::ddt(alphaPtr_(), precPKStar_[precI])
            +   fvm::Sp(dimensionedScalar("",dimless/dimTime,lambdas_[precI])*alphaPtr_(), precPKStar_[precI])
            -   isFuel_*initOneGroupFluxN_*dimensionedScalar("",dimPower,power_)*betas_[precI]/
            dimensionedScalar("",dimTime,promptGenerationTime_)
            +   fvm::div(phiPtr_(), precPKStar_[precI], "div(phi_,precStar_)")
            -   fvm::laplacian(diffCoeffPrecPtr_(), precPKStar_[precI])
        );

        scalar initRes = precEq.solve().max().initialResidual();
        neutroResidual = max(neutroResidual,initRes);

        //- Conver to precs per total volume
        precPK_[precI] = precPKStar_[precI] * alphaPtr_();
        precPK_[precI].correctBoundaryConditions();



        //- Calculate weighted precursors to be fed to power equation
        weightedPrec.set(precI,new scalar(fvc::domainIntegrate(
        initOneGroupFluxN_*precPK_[precI]).value()/domainIntegratedInitOneGroupFluxN_ ));

    }


    //- Calc power
    scalar sumPrecLamda = 0.0;
    forAll(precPK_,precI)
    {
        sumPrecLamda += weightedPrec[precI] * lambdas_[precI];
    }

    power_ = (powerOld_/dt+sumPrecLamda)/
    (
    1.0/dt - (totalReactivity_-beta_)/promptGenerationTime_
    );


    Info<<"Intermediate Power for the next iteration in the precursors equation[W] =  "<< power_ <<nl <<endl;
          
    neutroIter++;


    //- Some info on iteartions
    Info<< "neutroResidual (neutroIter: " << neutroIter << "): " << neutroResidual
    <<  endl;
}
while((neutroResidual > neutronIterationResidual)
      &&  (neutroIter < maxNeutronIterations)
    );

//- Re-scale fluxes, powerDensity
scalar fluxScaleFactor = power_/powerOld_;


for (int i = 0; i < energyGroups_; i++)
{
    fluxes_[i] *= fluxScaleFactor;
    fluxes_[i].correctBoundaryConditions();
}

oneGroupFlux_ *= fluxScaleFactor;
oneGroupFlux_.correctBoundaryConditions();

powerDensity_ *= fluxScaleFactor;
powerDensity_.correctBoundaryConditions();
secondaryPowerDensity_ *= fluxScaleFactor;
secondaryPowerDensity_.correctBoundaryConditions();

if(nuclearData_.found("updatePowerFields"))
{
    dictionary& subDict = nuclearData_.subDict("updatePowerFields");
    wordList regions = subDict.toc();
    forAll(regions, regI)
    {

        word fieldName = subDict.get<word>(regions[regI]);
        const dynamicFvMesh& regionMesh = mesh_.time().lookupObject<dynamicFvMesh>(regions[regI]);
        volScalarField& powerField = const_cast<volScalarField&>(regionMesh.lookupObject<volScalarField>(fieldName));
        powerField *= fluxScaleFactor;
    }
}



Info << endl << "pointKinetics: " << endl;
#include "pointKineticsInfo.H"
