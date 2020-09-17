//- Update reactivity coefficients
#include "computeFeedbackFieldValues.H"
#include "correctReactivity.H"

//- Update old time step values if time step has changed
if(timeIndex_ != mesh_.time().timeIndex())
{
    timeIndex_ = mesh_.time().timeIndex();
    powerOld_ = power_;
    //for (int i=0; i<delayedGroups_; i++)
    //{
    //    precursorPowersOld_[i] = precursorPowers_[i];
    //}
}

label n(delayedGroups_+1);
scalar dt(mesh_.time().deltaT().value());

//- Coefficient matrix, source, solution field
SquareMatrix<scalar> A(n, 0.0);
List<scalar> B(n, 0.0);
List<scalar> x(n, 0.0);

// calc precursors

PtrList<scalar> weightedPrec(delayedGroups_);
scalar domainIntegratedFluxSqr
(
    fvc::domainIntegrate(oneGroupFluxSqr).value()
);
volScalarField normaizedOneGroupFlux
(
    oneGroupFlux_/fvc::domainIntegrate(oneGroupFlux_)
);

    
for (int i = 0; i < 10; i++)
{
    bool initialize = false;
    if(initPrecursorsLiquidFuel_ &&  
        (mesh_.time().timeIndex() == mesh_.time().startTimeIndex()+1) && i==0)
    {
        initialize = true;
    }

    forAll(precPK_,precI)
    {
        if(initialize)
        {
            fvScalarMatrix precEq
            (
                fvm::Sp(dimensionedScalar("",dimless/dimTime,lambdas_[precI])*alphaPtr_(), precPKStar_[precI])
            -   normaizedOneGroupFlux*dimensionedScalar("",dimPower,power_)*betas_[precI]/
                dimensionedScalar("",dimTime,promptGenerationTime_)
            +   fvm::div(phiPtr_(), precPKStar_[precI], "div(phi_,precStar_)")
            -   fvm::laplacian(diffCoeffPrecPtr_(), precPKStar_[precI])                
            );
            precEq.solve();
            precPKStar_[precI].storePrevIter();  
        }
        else
        {
            fvScalarMatrix precEq
            (
                fvm::ddt(alphaPtr_(), precPKStar_[precI])
            +   fvm::Sp(dimensionedScalar("",dimless/dimTime,lambdas_[precI])*alphaPtr_(), precPKStar_[precI])
            -   normaizedOneGroupFlux*dimensionedScalar("",dimPower,power_)*betas_[precI]/
                dimensionedScalar("",dimTime,promptGenerationTime_)
            +   fvm::div(phiPtr_(), precPKStar_[precI], "div(phi_,precStar_)")
            -   fvm::laplacian(diffCoeffPrecPtr_(), precPKStar_[precI])                
            );
            precEq.solve();            
        }
        precPK_[precI] = precPKStar_[precI] * alphaPtr_(); // Precs per total volume
        precPK_[precI].correctBoundaryConditions();
        weightedPrec.set(precI,new scalar(fvc::domainIntegrate(
            oneGroupFluxSqr*precPK_[precI]).value()/domainIntegratedFluxSqr));
    }

    if(initialize)
    {
        precEquilibriumReactivity_ += beta_;
        forAll(precPK_,precI)
        {
            precEquilibriumReactivity_ -= weightedPrec[precI] * lambdas_[precI] * promptGenerationTime_ / powerOld_;
        }
    }


    //- First line of matrix has the coefficients of the discretized power 
    //  differential equation, the first element of the solution vector is the
    //  power
    scalar sumPrecLamda = 0.0;
    forAll(precPK_,precI)
    {
        sumPrecLamda += weightedPrec[precI] * lambdas_[precI];
    }
    power_ = (powerOld_/dt+sumPrecLamda)/
            (
                1.0/dt
                - (totalReactivity_-beta_)/promptGenerationTime_
            );

}
//- Re-scale fluxes, volFuelPower, precursors
scalar fluxScaleFactor = power_/powerOld_;

for (int i = 0; i < energyGroups_; i++)
{
    fluxes_[i] *= fluxScaleFactor;
    fluxes_[i].correctBoundaryConditions();
}
oneGroupFlux_ *= fluxScaleFactor;
oneGroupFlux_.correctBoundaryConditions();

volFuelPower_ *= fluxScaleFactor;
volFuelPower_.correctBoundaryConditions();

//- Only possible to re-scale precursors if they were initially present,
//  so I'm using a forAllIter rather than a for, as delayedGroups_ is not
//  necessarily representative of the length of the precursors_ PtrList
int i = 0;
forAllIter
(
    PtrList<volScalarField>,
    precursors_,
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
    reactorState_.set("precursorPowers", precursorPowers_);
    reactorState_.regIOobject::writeObject
    (
        IOstream::ASCII,
        IOstream::currentVersion,
        reactorState_.time().writeCompression(),
        true
    );
}

Info << endl << "pointKinetics: " << endl;
#include "pointKineticsInfo.H"
