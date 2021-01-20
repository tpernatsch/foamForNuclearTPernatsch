//- Get solution contro
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
#include "correctReactivity.H"

//- Update old time step values if time step has changed
if(timeIndex_ != mesh_.time().timeIndex())
{
    timeIndex_ = mesh_.time().timeIndex();
    powerOld_ = power_;
}

//label n(delayedGroups_+1);
scalar dt(mesh_.time().deltaT().value());

//- prepare fields necessary to weight precursors
PtrList<scalar> weightedPrec(delayedGroups_);
scalar domainIntegratedFluxSqr
(
    fvc::domainIntegrate(oneGroupFluxSqr).value()
);
volScalarField normaizedOneGroupFlux
(
    oneGroupFlux_/fvc::domainIntegrate(oneGroupFlux_)
);


//- Iterate precursors and power till sufficiently small prec residual
int neutroIter = 0;
scalar neutroResidual = 0.0;
do
{
    Info << "Neutron iteration no:  " << neutroIter << nl << endl;

    //- Check if you need to initialize precurosr distribution
    bool initialize = false;
    if(initPrecursorsLiquidFuel_ &&  
        (mesh_.time().timeIndex() == mesh_.time().startTimeIndex()+1) && neutroIter==0)
    {
        initialize = true;
    }

    //- Solve for precursors
    neutroResidual = 0.0;
    forAll(precPK_,precI)
    {
        //- Initialize precurosr distribution, if needed
        if(initialize)
        {
            fvMatrix<scalar> precEq
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

        //- Solve precursors
        fvMatrix<scalar> precEq
        (
            fvm::ddt(alphaPtr_(), precPKStar_[precI])
        +   fvm::Sp(dimensionedScalar("",dimless/dimTime,lambdas_[precI])*alphaPtr_(), precPKStar_[precI])
        -   normaizedOneGroupFlux*dimensionedScalar("",dimPower,power_)*betas_[precI]/
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
            oneGroupFluxSqr_*precPK_[precI]).value()/domainIntegratedFluxSqr));
    }

    //- If you have calculated a new precursor concentration at the initialization step,
    //  correct reactivity to start from a steady state
    if(initialize)
    {
        precEquilibriumReactivity_ += beta_;
        forAll(precPK_,precI)
        {
            precEquilibriumReactivity_ -= weightedPrec[precI] * lambdas_[precI] * promptGenerationTime_ / powerOld_;
        }
        #include "correctReactivity.H"
    }


    //- Calc power
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
    neutroIter++;

    //- Some info on iteartions 
    Info<< "neutroResidual (couplingIter: " << couplingIter 
        << ", neutroIter: " << neutroIter << "): " << neutroResidual 
        <<  endl;     
}
while
    ( 
        (neutroResidual > neutronIterationResidual) 
    &&  (neutroIter < maxNeutronIterations) 
    );

//- Re-scale fluxes, volFuelPower
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

Info << endl << "pointKinetics: " << endl;
#include "pointKineticsInfo.H"
