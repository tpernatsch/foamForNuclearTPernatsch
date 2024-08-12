/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "LimbackCreepModelLOCA.H"
#include "addToRunTimeSelectionTable.H"
#include "globalFieldLists.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * *//

namespace Foam
{
    defineTypeNameAndDebug(LimbackCreepModelLOCA, 0);
    addToRunTimeSelectionTable
    (
        creepModel, 
        LimbackCreepModelLOCA, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * ** * * * * //

Foam::LimbackCreepModelLOCA::LimbackCreepModelLOCA
(
    const fvMesh& mesh,
    const dictionary& lawDict
)
:
    creepModel(mesh, lawDict),
    fastFluenceName_(lawDict.lookupOrDefault<word>("fastFluenceName", "fastFluence" )),
    fastFluxName_(lawDict.lookupOrDefault<word>("fastFluxName", "fastFlux" )),
    fastFluence_(0),
    fastFlux_(0),
    betaFraction_(0),
    epsilonCreepIrr_(createOrLookup<symmTensor>(mesh, "epsilonCreepIrr")),
    epsilonCreepTh_(createOrLookup<symmTensor>(mesh, "epsilonCreepTh")),
    epsilonCreepPrim_(createOrLookup<symmTensor>(mesh, "epsilonCreepPrim")),
    epsilonCreepIrrEq_(createOrLookup<scalar>(mesh, "epsilonCreepIrrEq")),
    DepsilonCreepIrrEq_(createOrLookup<scalar>(mesh, "DepsilonCreepIrrEq")),
    epsilonCreepThEq_(createOrLookup<scalar>(mesh, "epsilonCreepThEq")),
    DepsilonCreepThEq_(createOrLookup<scalar>(mesh, "DepsilonCreepThEq")),
    epsilonCreepPrimEq_(createOrLookup<scalar>(mesh, "epsilonCreepPrimEq")),
    DepsilonCreepPrimEq_(createOrLookup<scalar>(mesh, "DepsilonCreepPrimEq")),
    timeIndex_(0),
    relax_(lawDict.lookupOrDefault<scalar>("relax", 1.0 )),
    irradiationCreep_(lawDict.lookupOrDefault<bool>("irradiationCreep", true)),
    primaryCreep_(lawDict.lookupOrDefault<bool>("primaryCreep", true)),
    NewtonRaphsonMethod_(lawDict.lookupOrDefault<bool>("NewtonRaphsonMethod", false)),
    NewtonRaphsonTolerance_(NewtonRaphsonMethod_ ? 
        lawDict.lookupOrDefault<scalar>("NewtonRaphsonTolerance", 1e-6) : 1.0),
    C0(3.557e-24),
    C1(0.85),
    C2(1.0)
{
    epsilonCreepIrr_.oldTime();
    epsilonCreepTh_.oldTime();
    epsilonCreepPrim_.oldTime();
    epsilonCreepIrrEq_.oldTime();
    DepsilonCreepIrrEq_.oldTime();
    epsilonCreepThEq_.oldTime();
    DepsilonCreepThEq_.oldTime();
    epsilonCreepPrimEq_.oldTime();
    DepsilonCreepPrimEq_.oldTime();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::LimbackCreepModelLOCA::~LimbackCreepModelLOCA()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::LimbackCreepModelLOCA::correctCreep
(
    volSymmTensorField& epsilonEl,
    const labelList& addr
)
{
    if(fastFluence_ == nullptr)
    {        
        if(mesh_.foundObject<volScalarField>(fastFluenceName_))
        {
            fastFluence_ = &mesh_.lookupObject<volScalarField>(fastFluenceName_);
        }
        else
        {
            FatalErrorInFunction
                << nl
                << "    Fast fluence field \"" << fastFluenceName_ << "\" requested by " 
                << "LimbackCreepModelLOCA but not found in registry." << nl
                << abort(FatalError);
        }

        if(mesh_.foundObject<volScalarField>(fastFluxName_))
        {
            fastFlux_ = &mesh_.lookupObject<volScalarField>(fastFluxName_);
        }
        else
        {
            FatalErrorInFunction
                << nl
                << "    Fast flux field \"" << fastFluxName_ << "\" requested by " 
                << "LimbackCreepModelLOCA but not found in registry." << nl
                << abort(FatalError);
        }
    }  

    if(betaFraction_ == nullptr)
    {        
        if(mesh_.foundObject<volScalarField>("betaFraction"))
        {
            betaFraction_ = &mesh_.lookupObject<volScalarField>("betaFraction");
        }
        else
        {
            FatalErrorInFunction
                << nl
                << "    Beta fraction field \"betaFraction\" requested by " 
                << "LimbackCreepModelLOCA but not found in registry." << nl
                << abort(FatalError);
        }
    } 

    const volScalarField& mu(mesh_.lookupObject<volScalarField>("mu"));
    const volScalarField& lambda(mesh_.lookupObject<volScalarField>("lambda"));
    const volScalarField& Temp(mesh_.lookupObject<volScalarField>("T"));

    // Initialize tolerance for inner convergence loop
    scalar tol(1.0);

    // Irradiation time step [hr]
    scalar deltaT(0);
    
    if(mesh_.time().timeIndex() > 0)
    {
        deltaT = mesh_.time().deltaTValue()/3600.0;
        
        // Initialize creep thermal rate from previous iterations
        volScalarField creepThermalRate(DepsilonCreepThEq_/deltaT);

        const scalarField& epsilonCreepIrrEqOld_ = 
        epsilonCreepIrrEq_.oldTime().internalField();
        const scalarField& epsilonCreepThEqOld_ = 
        epsilonCreepThEq_.oldTime().internalField();
        const scalarField& epsilonCreepPrimEqOld_ = 
        epsilonCreepPrimEq_.oldTime().internalField();

        const scalarField& fastFluenceI_ = fastFluence_->internalField();
        const scalarField& fastFluxI_ = fastFlux_->internalField();

        const scalarField& epsilonCreepEqOld_ = 
        epsilonCreepEq_.oldTime().internalField();

        const symmTensorField& epsilonCreepIrrOldI = 
        epsilonCreepIrr_.oldTime().internalField();
        const symmTensorField& epsilonCreepThOldI = 
        epsilonCreepTh_.oldTime().internalField();
        const symmTensorField& epsilonCreepPrimOldI = 
        epsilonCreepPrim_.oldTime().internalField();
        const symmTensorField& epsilonCreepOldI = 
        epsilonCreep_.oldTime().internalField();
    
        // Adjust elastic strain tensor by subtracting old creep values    
        forAll(addr, addrI)
        {
            const label cellI = addr[addrI];

            epsilonEl[cellI] -= epsilonCreepOldI[cellI];

            const cell& c = mesh_.cells()[cellI];     

            forAll(c, faceI)
            {
                const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);

                {
                    if (patchID > -1 and epsilonEl.boundaryField()[patchID].size())
                    {
                        const label faceID = 
                        mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);
                        
                        symmTensorField& epsilonElP = 
                        epsilonEl.boundaryFieldRef()[patchID];
                        const symmTensorField& epsilonCreepOldP = 
                        epsilonCreep_.oldTime().boundaryField()[patchID];

                        epsilonElP[faceID] -= epsilonCreepOldP[faceID];
                    }
                }   
            }
        }    

        // If in first iteration of new time step, set creep strain increment
        // to 0  
        if(mesh_.time().timeIndex() > timeIndex_)
        { 
            forAll(addr, addrI)
            {
                const label cellI = addr[addrI];
                DepsilonCreepEq_.ref()[cellI] *= 0;
                DepsilonCreepIrrEq_.ref()[cellI] *= 0;
                DepsilonCreepThEq_.ref()[cellI] *= 0;
                DepsilonCreepPrimEq_.ref()[cellI] *= 0;
                
                creepThermalRate[cellI] *= 0;

                const cell& c = mesh_.cells()[cellI];  

                forAll(c, faceI)
                {
                    const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);

                    {
                        if (patchID > -1 and epsilonEl.boundaryField()[patchID].size())
                        {
                            const label faceID = 
                            mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);

                            scalarField& DepsilonCreepEqP = 
                            DepsilonCreepEq_.boundaryFieldRef()[patchID];
                            scalarField& DepsilonCreepIrrEqP = 
                            DepsilonCreepIrrEq_.boundaryFieldRef()[patchID];
                            scalarField& DepsilonCreepThP = 
                            DepsilonCreepThEq_.boundaryFieldRef()[patchID];
                            scalarField& DepsilonCreepPrimEqP = 
                            DepsilonCreepPrimEq_.boundaryFieldRef()[patchID];

                            scalarField& creepThermalRateP = 
                            creepThermalRate.boundaryFieldRef()[patchID];

                            DepsilonCreepEqP[faceID] *= 0;
                            DepsilonCreepIrrEqP[faceID] *= 0;
                            DepsilonCreepThP[faceID] *= 0;
                            DepsilonCreepPrimEqP[faceID] *= 0;

                            creepThermalRateP[faceID] *= 0;
                        }
                    }
                }
            }

            timeIndex_ = mesh_.time().timeIndex();
        }

        // Do-while loop to correct creep strain
        do
        {
            scalar maxRate(0);
            scalar maxFlux(0);
            scalar maxTemp(0);
            scalar maxStress(0);

            scalar creepIrrRate(0);
            scalar creepIrrRatePrime(0);

            scalar creepThRate(0);
            scalar creepThRatePrime(0);

            scalar creepPrimRate(0);
            scalar creepPrimRatePrime(0);

            forAll(addr, addrI)
            {
                const label cellI = addr[addrI];

                // Current stress (without latest creep strain increment)
                const symmTensor sigmaI = 
                2*mu[cellI]*epsilonEl[cellI] + lambda[cellI]*tr(epsilonEl[cellI])*I;

                // Calculate the von Mises stress [Pa] (considering latest creep
                //  increment)
                const scalar sigmaEff = 
                max(min
                (
                    sqrt((3.0/2.0)
                    *magSqr(dev(sigmaI))) 
                    - 3*mu[cellI]*DepsilonCreepEq_.ref()[cellI] 
                    , 3e10
                ),
                1e5);
                
                maxStress = max(sigmaEff, maxStress);

                const scalar DepsilonCreepEqPrev = DepsilonCreepEq_.ref()[cellI];

                // Reference to temperature field [K]
                const scalar T = Temp.internalField()[cellI];

                maxTemp = max(T, maxTemp);

                /****** Calculate irradiation creep rate [1/hr] ********/

                if (irradiationCreep_)
                {
                    scalar phi = fastFluxI_[cellI]; 
                    maxFlux = max(phi, maxFlux);

                    creepIrrRate = calcIrradiationCreepRate(phi, sigmaEff);

                    if (NewtonRaphsonMethod_)
                    {
                        creepIrrRatePrime = calcIrradiationCreepRatePrime(phi, sigmaEff, 
                            creepIrrRate);
                    }

                    // Update irradiation creep equivalent
                    DepsilonCreepIrrEq_.ref()[cellI] =  
                    (1 - relax_)*DepsilonCreepIrrEq_.ref()[cellI] 
                    + relax_*creepIrrRate*deltaT;

                    epsilonCreepIrrEq_.ref()[cellI] = 
                    epsilonCreepIrrEqOld_[cellI] + DepsilonCreepIrrEq_.ref()[cellI];
                }

                /****** Calculate secondary thermal creep rate [1/hr] ********/

                // Access beta fraction of cellI
                const scalar betaFrac = betaFraction_->internalField()[cellI];

                if( T <= 700 )
                {
                    // Compute thermal creep in normal operating conditions:
                    creepThRate = calcNormalThermalCreepRate(
                        fastFluenceI_[cellI], T, sigmaEff);

                    if (NewtonRaphsonMethod_)
                    {
                        // Compute derivative of thermal creep in normal operating conditions:
                        creepThRatePrime = calcNormalThermalCreepRatePrime(
                            fastFluenceI_[cellI], T, sigmaEff, creepThRate);
                    }
                }
                else if( T > 700 and T < 900)
                {
                    // Compute thermal creep in normal operating conditions:
                    scalar cTh700 = calcNormalThermalCreepRate(
                        fastFluenceI_[cellI], 700, sigmaEff);

                    // Compute thermal creep in LOCA conditions (Erbacher):
                    scalar cTh900 = calcLOCAThermalCreepRate(
                        900, betaFrac, sigmaEff, creepThermalRate[cellI]);

                    // Interpolation:
                    creepThRate   = cTh700 + (
                        (T-700)/(900-700) * (cTh900-cTh700));

                    if (NewtonRaphsonMethod_)
                    {
                        // Compute derivative of thermal creep in normal operating conditions:
                        scalar cThPrime700 = calcNormalThermalCreepRatePrime(
                            fastFluenceI_[cellI], 700, sigmaEff, creepThRate);

                        // Compute derivative of thermal creep in LOCA conditions (Erbacher):
                        scalar cThPrime900 = calcLOCAThermalCreepRatePrime(
                            900, betaFrac, sigmaEff, cTh900);

                        creepThRatePrime = cThPrime700 + (
                            (T-700)/(900-700) * (cThPrime900-cThPrime700));
                    }
                }
                else if( T >= 900)
                {
                    // Compute thermal creep in LOCA conditions (Erbacher):
                    creepThRate = calcLOCAThermalCreepRate(
                        T, betaFrac, sigmaEff, creepThermalRate[cellI]);

                    if (NewtonRaphsonMethod_)
                    {
                        // Compute derivative of thermal creep in LOCA conditions (Erbacher):
                        creepThRatePrime = calcLOCAThermalCreepRatePrime(
                            T, betaFrac, sigmaEff, creepThRate);
                    }
                }

                creepThermalRate[cellI] = creepThRate;

                // Update thermal creep equivalent
                DepsilonCreepThEq_.ref()[cellI] = 
                (1 - relax_)*DepsilonCreepThEq_.ref()[cellI] + 
                relax_*creepThRate*deltaT;

                epsilonCreepThEq_.ref()[cellI] = 
                epsilonCreepThEqOld_[cellI] + DepsilonCreepThEq_.ref()[cellI];

                /****** Calculate primary thermal creep rate ********/

                if (primaryCreep_)
                {                    
                    creepPrimRate = calcPrimaryCreepRate(creepIrrRate, creepThRate, 
                        epsilonCreepPrimEqOld_[cellI]);
                    
                    if (NewtonRaphsonMethod_)
                    {
                        creepPrimRatePrime = calcPrimaryCreepRatePrime
                        (
                            creepIrrRate, 
                            creepThRate, 
                            creepIrrRatePrime, 
                            creepThRatePrime,
                            epsilonCreepPrimEq_.ref()[cellI], 
                            epsilonCreepPrimEqOld_[cellI]
                        );
                    }

                    DepsilonCreepPrimEq_.ref()[cellI] = 
                    (1 - relax_)*DepsilonCreepPrimEq_.ref()[cellI] + 
                    relax_*creepPrimRate*deltaT;

                    epsilonCreepPrimEq_.ref()[cellI] = epsilonCreepPrimEqOld_[cellI] + 
                    DepsilonCreepPrimEq_.ref()[cellI];
                }

                /******* Total creep rate ********/ 

                const scalar fTrial
                (
                    creepIrrRate 
                  + creepThRate 
                  + creepPrimRate 
                );

                if (NewtonRaphsonMethod_)
                {
                    const scalar fTrialPrime
                    (
                        creepIrrRatePrime
                      + creepThRatePrime 
                      + creepPrimRatePrime 
                    );

                    DepsilonCreepEq_.ref()[cellI] = 
                    (1 - relax_)*DepsilonCreepEq_.ref()[cellI] 
                    + relax_*
                    (
                        DepsilonCreepEq_.ref()[cellI] - 
                        (
                            DepsilonCreepEq_.ref()[cellI] - deltaT*fTrial
                        )/(1 + 3*mu[cellI]/1e6*fTrialPrime*deltaT)
                    );   
                }
                else
                {
                    DepsilonCreepEq_.ref()[cellI] = 
                    (1 - relax_)*DepsilonCreepEq_.ref()[cellI] + relax_*fTrial*deltaT;
                }                             

                epsilonCreepEq_.ref()[cellI] = 
                epsilonCreepEqOld_[cellI] + DepsilonCreepEq_.ref()[cellI];
                
                maxRate = max(fTrial, maxRate);

                /******* Update 3D creep fields ********/ 

                epsilonCreepIrr_.ref()[cellI] = 
                epsilonCreepIrrOldI[cellI] + 1.5*DepsilonCreepIrrEq_.ref()[cellI]*dev(sigmaI 
                    - 2*mu[cellI]*DepsilonCreep_.ref()[cellI])/max(sigmaEff, SMALL);

                epsilonCreepTh_.ref()[cellI] = 
                epsilonCreepThOldI[cellI] + 1.5*DepsilonCreepThEq_.ref()[cellI]*dev(sigmaI  
                    - 2*mu[cellI]*DepsilonCreep_.ref()[cellI])/max(sigmaEff, SMALL); 

                epsilonCreepPrim_.ref()[cellI] = 
                epsilonCreepPrimOldI[cellI] + 1.5*DepsilonCreepPrimEq_.ref()[cellI]*dev(sigmaI 
                    - 2*mu[cellI]*DepsilonCreep_.ref()[cellI])/max(sigmaEff, SMALL);

                DepsilonCreep_.ref()[cellI] = 
                (1 - relax_)*DepsilonCreep_.ref()[cellI] 
                + relax_*1.5*DepsilonCreepEq_.ref()[cellI]*dev(sigmaI 
                    - 2*mu[cellI]*DepsilonCreep_.ref()[cellI])/max(sigmaEff, SMALL); 

                epsilonCreep_.ref()[cellI] = 
                epsilonCreepOldI[cellI] + DepsilonCreep_.ref()[cellI]; 

                const scalar tolTr = min(1, mag( (DepsilonCreepEq_.ref()[cellI] - DepsilonCreepEqPrev)
                    / (max(DepsilonCreepEqPrev , SMALL))));

                tol = min(tol, tolTr );

                const cell& c = mesh_.cells()[cellI];  

                // Loop over all faces of current cell
                forAll(c, faceI)
                {
                    const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);
                    {
                        if (patchID > -1 and epsilonEl.boundaryField()[patchID].size())
                        {
                            const label faceID = 
                            mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);

                            const scalarField& Tp =
                            Temp.boundaryField()[patchID];
                            const scalarField& muP = 
                            mu.boundaryField()[patchID];
                            const scalarField& lambdaP = 
                            lambda.boundaryField()[patchID];
                            const scalarField& fastFluenceP = 
                            fastFluence_->boundaryField()[patchID];
                            const scalarField& fastFluxP = 
                            fastFlux_->boundaryField()[patchID];
                            const symmTensorField& epsilonElP = 
                            epsilonEl.boundaryField()[patchID];

                            symmTensorField& epsilonCreepP = 
                            epsilonCreep_.boundaryFieldRef()[patchID];
                            symmTensorField& DepsilonCreepP = 
                            DepsilonCreep_.boundaryFieldRef()[patchID];
                            symmTensorField& epsilonCreepIrrP = 
                            epsilonCreepIrr_.boundaryFieldRef()[patchID];
                            symmTensorField& epsilonCreepThP = 
                            epsilonCreepTh_.boundaryFieldRef()[patchID];
                            symmTensorField& epsilonCreepPrimP = 
                            epsilonCreepPrim_.boundaryFieldRef()[patchID];

                            scalarField& epsilonCreepEqP = 
                            epsilonCreepEq_.boundaryFieldRef()[patchID];
                            scalarField& epsilonCreepIrrEqP = 
                            epsilonCreepIrrEq_.boundaryFieldRef()[patchID];
                            scalarField& epsilonCreepThEqP = 
                            epsilonCreepThEq_.boundaryFieldRef()[patchID];
                            scalarField& epsilonCreepPrimEqP = 
                            epsilonCreepPrimEq_.boundaryFieldRef()[patchID];

                            const scalarField& epsilonCreepIrrEqOldP = 
                            epsilonCreepIrrEq_.oldTime().boundaryField()[patchID];
                            const scalarField& epsilonCreepThEqOldP = 
                            epsilonCreepThEq_.oldTime().boundaryField()[patchID];
                            const scalarField& epsilonCreepPrimEqOldP = 
                            epsilonCreepPrimEq_.oldTime().boundaryField()[patchID];
                            const scalarField& epsilonCreepEqOldP = 
                            epsilonCreepEq_.oldTime().boundaryField()[patchID];

                            const symmTensorField& epsilonCreepIrrOldP = 
                            epsilonCreepIrr_.oldTime().boundaryField()[patchID];
                            const symmTensorField& epsilonCreepThOldP = 
                            epsilonCreepTh_.oldTime().boundaryField()[patchID];
                            const symmTensorField& epsilonCreepPrimOldP = 
                            epsilonCreepPrim_.oldTime().boundaryField()[patchID];
                            const symmTensorField& epsilonCreepOldP = 
                            epsilonCreep_.oldTime().boundaryField()[patchID];

                            scalarField& DepsilonCreepIrrEqP = 
                            DepsilonCreepIrrEq_.boundaryFieldRef()[patchID];
                            scalarField& DepsilonCreepThEqP = 
                            DepsilonCreepThEq_.boundaryFieldRef()[patchID];
                            scalarField& DepsilonCreepPrimEqP = 
                            DepsilonCreepPrimEq_.boundaryFieldRef()[patchID];
                            scalarField& DepsilonCreepEqP = 
                            DepsilonCreepEq_.boundaryFieldRef()[patchID];

                            scalarField& creepThermalRateP = 
                            creepThermalRate.boundaryFieldRef()[patchID];

                            // Calculate the von Mises stress [Pa]
                            const symmTensor sigmaPI = 
                            2*muP[faceID]*epsilonElP[faceID] + 
                            lambdaP[faceID]*tr(epsilonElP[faceID])*I;

                            const scalar sigmaEff = 
                            max
                            (
                                min
                                (
                                    sqrt((3.0/2.0)*magSqr(dev(sigmaPI))) 
                                    - 3*muP[faceID]*DepsilonCreepEqP[faceID], 
                                    3e10
                                )
                                , 1e5
                            );

                            const scalar DepsilonCreepEqPPrev = 
                            DepsilonCreepEqP[faceID];

                            // Reference to temperature field [K]
                            const scalar T = Tp[faceID];                    

                            /****** Calculate irradiation creep rate [1/hr]  *********/   

                            if (irradiationCreep_)
                            {   
                                // Fast flux 
                                scalar phi = fastFluxP[faceID];

                                creepIrrRate = calcIrradiationCreepRate(phi, sigmaEff);

                                if (NewtonRaphsonMethod_)
                                {
                                    creepIrrRatePrime = calcIrradiationCreepRatePrime(phi, sigmaEff, 
                                        creepIrrRate);
                                }
                                
                                DepsilonCreepIrrEqP[faceID] =  
                                (1 - relax_)*DepsilonCreepIrrEqP[faceID] 
                                + relax_*creepIrrRate*deltaT;

                                epsilonCreepIrrEqP[faceID] = 
                                DepsilonCreepIrrEqP[faceID] + epsilonCreepIrrEqOldP[faceID];
                            }

                        /****** Calculate secondary thermal creep rate [1/hr] ********/

                            // Access beta fraction of cellI
                            const scalar betaFrac = betaFraction_->boundaryField()[patchID][faceID];

                            if ( T <= 700 )
                            {
                                // Compute thermal creep in normal operating conditions:
                                creepThRate = calcNormalThermalCreepRate(
                                    fastFluenceP[faceID], T, sigmaEff);

                                if (NewtonRaphsonMethod_)
                                {
                                    // Compute derivative of thermal creep in normal operating conditions:
                                    creepThRatePrime = calcNormalThermalCreepRatePrime(
                                        fastFluenceP[faceID], T, sigmaEff, creepThRate);
                                }
                            }
                            else if( T > 700 and T < 900)
                            {
                                // Compute thermal creep in normal operating conditions:
                                scalar cTh700 = calcNormalThermalCreepRate(
                                    fastFluenceP[faceID], 700, sigmaEff);

                                // Compute thermal creep in LOCA conditions (Erbacher):
                                scalar cTh900 = calcLOCAThermalCreepRate(
                                    900, betaFrac, sigmaEff, creepThermalRateP[faceID]);

                                // Interpolation:
                                creepThRate   = cTh700 + (
                                    (T-700)/(900-700) * (cTh900-cTh700));

                                if (NewtonRaphsonMethod_)
                                {
                                    // Compute derivative of thermal creep in normal operating conditions:
                                    scalar cThPrime700 = calcNormalThermalCreepRatePrime(
                                        fastFluenceP[faceID], 700, sigmaEff, cTh700);

                                    // Compute derivative of thermal creep in LOCA conditions (Erbacher):
                                    scalar cThPrime900 = calcLOCAThermalCreepRatePrime(
                                        900, betaFrac, sigmaEff, cTh900);

                                    creepThRatePrime = cThPrime700 + (
                                        (T-700)/(900-700) * (cThPrime900-cThPrime700));
                                }
                            }
                            else if ( T >= 900 )
                            {
                                // Compute thermal creep in LOCA conditions (Erbacher):
                                creepThRate = calcLOCAThermalCreepRate(
                                    T, betaFrac, sigmaEff, creepThermalRateP[faceID]);

                                if (NewtonRaphsonMethod_)
                                {
                                    // Compute derivative of thermal creep in LOCA conditions (Erbacher):
                                    creepThRatePrime = calcLOCAThermalCreepRatePrime(
                                        T, betaFrac, sigmaEff, creepThRate);
                                }
                            }

                            creepThermalRateP[faceID] = creepThRate;

                            // Update thermal creep equivalent
                            DepsilonCreepThEqP[faceID] = 
                            (1 - relax_)*DepsilonCreepThEqP[faceID] + 
                            relax_*creepThRate*deltaT;               

                            epsilonCreepThEqP[faceID] = 
                            DepsilonCreepThEqP[faceID] + epsilonCreepThEqOldP[faceID] ;

                            /****** Calculate primary thermal creep rate ********/

                            if (primaryCreep_)
                            {                
                                creepPrimRate = calcPrimaryCreepRate(creepIrrRate, creepThRate, 
                                    epsilonCreepPrimEqOldP[faceID]);

                                if (NewtonRaphsonMethod_)
                                {
                                    creepPrimRatePrime = calcPrimaryCreepRatePrime
                                    (
                                        creepIrrRate, 
                                        creepThRate, 
                                        creepIrrRatePrime, 
                                        creepThRatePrime,
                                        epsilonCreepPrimEqP[faceID], 
                                        epsilonCreepPrimEqOldP[faceID]
                                    ); 
                                }

                                DepsilonCreepPrimEqP[faceID] = 
                                (1 - relax_)*DepsilonCreepPrimEqP[faceID] + 
                                relax_*creepPrimRate*deltaT;

                                epsilonCreepPrimEqP[faceID] = 
                                epsilonCreepPrimEqOldP[faceID] + 
                                DepsilonCreepPrimEqP[faceID];
                            }

                            /******* Total creep rate ********/

                            const scalar fTrial
                            (
                                creepIrrRate 
                              + creepThRate
                              + creepPrimRate 
                            );
                
                            maxRate = max(fTrial, maxRate);

                            if (NewtonRaphsonMethod_)
                            {
                                const scalar fTrialPrime
                                (
                                    creepIrrRatePrime
                                  + creepThRatePrime 
                                  + creepPrimRatePrime 
                                );

                                DepsilonCreepEqP[faceID] = 
                                (1 - relax_)*DepsilonCreepEqP[faceID] 
                                + relax_*
                                (
                                    DepsilonCreepEqP[faceID] - 
                                    (
                                        DepsilonCreepEqP[faceID] - deltaT*fTrial
                                    )/(1 + 3*muP[faceID]/1e6*fTrialPrime*deltaT)
                                );
                            }
                            else
                            {
                                DepsilonCreepEqP[faceID] = 
                                (1-relax_)*DepsilonCreepEqP[faceID] + relax_*fTrial*deltaT;                              
                            }

                            epsilonCreepEqP[faceID] = 
                            DepsilonCreepEqP[faceID] + epsilonCreepEqOldP[faceID];

                            /******* Update 3D creep fields ********/ 
                            
                            epsilonCreepIrrP[faceID] = 
                            epsilonCreepIrrOldP[faceID] + 1.5*DepsilonCreepIrrEqP[faceID]*dev(
                                sigmaPI - 2*muP[faceID]*DepsilonCreepP[faceID])/max(sigmaEff, SMALL);
                            
                            epsilonCreepThP[faceID] = 
                            epsilonCreepThOldP[faceID] + 1.5*DepsilonCreepThEqP[faceID]*dev(
                                sigmaPI - 2*muP[faceID]*DepsilonCreepP[faceID])/max(sigmaEff, SMALL);

                            epsilonCreepPrimP[faceID] = 
                            epsilonCreepPrimOldP[faceID] + 1.5*DepsilonCreepPrimEqP[faceID]*dev(
                                sigmaPI - 2*muP[faceID]*DepsilonCreepP[faceID])/max(sigmaEff, SMALL);

                            DepsilonCreepP[faceID] = 
                            (1 - relax_)*DepsilonCreepP[faceID] + relax_*1.5*DepsilonCreepEqP[faceID]*dev(
                                sigmaPI - 2*muP[faceID]*DepsilonCreepP[faceID])/max(sigmaEff, SMALL) ;
                            
                            epsilonCreepP[faceID] = 
                            epsilonCreepOldP[faceID] + DepsilonCreepP[faceID];
                            
                            const scalar tolTr = 
                            min
                            (
                                1, 
                                mag
                                ( 
                                    (DepsilonCreepEqP[faceID] - DepsilonCreepEqPPrev)/ 
                                    (max(DepsilonCreepEqPPrev , SMALL))
                                )
                            );

                            tol = min(tol, tolTr );
                        }
                    }
                }
            }
            
            // Info << "Limback, maxTemp: " << maxTemp << ", maxFlux: " << maxFlux << ", maxStress: " << maxStress << ", maxRate (1/hr) " << maxRate << endl;

        }while(tol > NewtonRaphsonTolerance_);
    }    
    
    forAll(addr, addrI)
    {
        const label cellI = addr[addrI];

        epsilonEl[cellI] -= DepsilonCreep_.ref()[cellI];     
        
        const cell& c = mesh_.cells()[cellI];         

        forAll(c, faceI)
        {
            const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);   

            {
                if (patchID > -1 and epsilonEl.boundaryField()[patchID].size())
                {
                    const label faceID = 
                    mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);
                    
                    symmTensorField& epsilonElP = 
                    epsilonEl.boundaryFieldRef()[patchID];
                    const symmTensorField& DepsilonCreepP = 
                    DepsilonCreep_.boundaryField()[patchID];

                    epsilonElP[faceID] -= DepsilonCreepP[faceID];
                }
            }   
        }
    }
}


Foam::scalar Foam::LimbackCreepModelLOCA::calcIrradiationCreepRate
(
    const scalar& flux,
    const scalar& sigmaEff
)
{
    const scalar C0(3.557e-24);
    const scalar C1(0.85);
    const scalar C2(1.0);

    return C0*pow(flux*1e4, C1)*pow(sigmaEff/1e6, C2);
}


Foam::scalar Foam::LimbackCreepModelLOCA::calcIrradiationCreepRatePrime
(
    const scalar& flux,
    const scalar& sigmaEff,
    const scalar& creepIrrRate
)
{
    const scalar C2(1.0);

    return C2*creepIrrRate/(sigmaEff/1e6);
}


Foam::scalar Foam::LimbackCreepModelLOCA::calcNormalThermalCreepRate
(
    const scalar& fastFluence,
    const scalar& T,
    const scalar& sigmaEff
)
{
    // Gas Constant [ J/mol/K]
    const scalar R = 8.314;
    // Material constant [K/Mpa/hr]
    const scalar A = 1.08e9;
    // Activation energy per mol [J/mol]
    const scalar Q = 201e3;
    // Primary creep exponent 
    const scalar n = 2;
    // Dimensionless constants
    const scalar A1 = 0.56;
    const scalar a = 650;
    const scalar A3 = 1.3;
    // Material constant [n/cm2]^-A3
    const scalar A2 = 1.4e-27;

    const scalar E = 1.148e5 - 59.9*T;

    const scalar ai = a*(1 - A1*(1 - exp(-A2*pow(fastFluence,A3))));

    scalar creepThRate = 
    A*E/T*pow(Foam::sinh(ai*sigmaEff/1e6/E), n)*exp(-Q/R/T);

    return creepThRate;
}


Foam::scalar Foam::LimbackCreepModelLOCA::calcNormalThermalCreepRatePrime
(
    const scalar& fastFluence,
    const scalar& T,
    const scalar& sigmaEff,
    const scalar& creepThRate
)
{
    // Primary creep exponent 
    const scalar n = 2;
    // Dimensionless constants
    const scalar A1 = 0.56;
    const scalar a = 650;
    const scalar A3 = 1.3;
    // Material constant [n/cm2]^-A3
    const scalar A2 = 1.4e-27;

    const scalar E = 1.148e5 - 59.9*T;

    const scalar ai = a*(1 - A1*(1 - exp(-A2*pow(fastFluence,A3))));

    return n*ai*creepThRate/E*(1.0/Foam::tanh(ai*sigmaEff/1e6/E));
}


Foam::scalar Foam::LimbackCreepModelLOCA::calcPrimaryCreepRate
(
    const scalar& creepIrrRate,
    const scalar& creepThRate,
    const scalar& epsilonCreepPrimEqOld
)
{
    // Dimensionless constants
    const scalar C(52);
    const scalar b(0.109);
    const scalar d(-2.05);
    //Maerial constant [hr^b]
    const scalar B(0.0216);
    //Maerial constant [hr]
    const scalar D(35500);

    // DeltaT
    const scalar deltaT = mesh_.time().deltaTValue()/3600.0;

    // Define secondary creep rate
    const scalar creepSecRate(creepIrrRate + creepThRate);

    if(creepSecRate > VSMALL)
    {
        // Saturated primary creep
        const scalar epsilonPrim_sat = 
        B*pow(creepSecRate, b)*pow(2 - Foam::tanh(D*creepSecRate), d);

        if(epsilonCreepPrimEqOld < (epsilonPrim_sat-SMALL))
        {
            // Time constant
            const scalar tau = 
            1/C/C*pow
            (
                Foam::log
                (
                    1 - epsilonCreepPrimEqOld/epsilonPrim_sat
                )
                ,2
            )
            *1/creepSecRate;

            // Compute and return primary creep rate
            const scalar creepPrimRate = 1/deltaT*(epsilonPrim_sat*(1 - exp(
                -C*sqrt(creepSecRate*(tau + deltaT)))) - epsilonCreepPrimEqOld);

            return creepPrimRate;
        }
    }
    
    return 0;
}


Foam::scalar Foam::LimbackCreepModelLOCA::calcPrimaryCreepRatePrime
(
    const scalar& creepIrrRate,
    const scalar& creepThRate,
    const scalar& creepIrrRatePrime,
    const scalar& creepThRatePrime,
    const scalar& epsilonCreepPrimEq,
    const scalar& epsilonCreepPrimEqOld
)
{
    // Dimensionless constants
    const scalar C(52);
    const scalar b(0.109);
    const scalar d(-2.05);
    //Maerial constant [hr^b]
    const scalar B(0.0216);
    //Maerial constant [hr]
    const scalar D(35500);

    // DeltaT
    const scalar deltaT = mesh_.time().deltaTValue()/3600.0;

    // Define secondary creep rate
    const scalar creepSecRate(creepIrrRate + creepThRate);
    const scalar creepSecRatePrime(creepIrrRatePrime + creepThRatePrime);

    if(creepSecRate > VSMALL)
    {
        // Saturated primary creep
        const scalar epsilonPrim_sat = 
        B*pow(creepSecRate, b)*pow(2 - Foam::tanh(D*creepSecRate), d);

        if(epsilonCreepPrimEqOld < (epsilonPrim_sat-SMALL))
        {
            const scalar alpha = 
            max(1 - epsilonCreepPrimEq/epsilonPrim_sat, SMALL);
            const scalar alphaOld = 
            1 - epsilonCreepPrimEqOld/epsilonPrim_sat;

            scalar creepPrime = (creepSecRatePrime)*
            (
                1/deltaT*
                (
                    B/creepSecRate - 
                    (
                        d*D*pow(Foam::cosh(D*creepSecRate), -2.0)
                    )/(2 - Foam::tanh(D*creepSecRate))
                )*
                (
                    epsilonCreepPrimEq - 
                    alpha*Foam::log(alphaOld)*epsilonCreepPrimEqOld/min(alphaOld*Foam::log(alpha), -VSMALL)
                ) -
                epsilonPrim_sat*C*C*alpha/min(Foam::log(alpha), -VSMALL)
            );

            return creepPrime;
        }
    }

    return 0;
}


Foam::scalar Foam::LimbackCreepModelLOCA::calcLOCAThermalCreepRate
(
    const scalar& T,
    const scalar& betaFrac,
    const scalar& sigmaEff,
    scalar creepThRateOld
)
{
    // TODO (EB) : this model makes use of the betaFraction field. Currently 
    // OFFBEAT creates the field (and initializes it to zero) even if the model
    // is explicitly set to "none" (coeherently to what swellingModel does). 
    // Now this is not a problem, as in this case the correlations will be 
    // evaluated as the material was pure alpha. In the future, if it is
    // decided that the field betaFraction is created only when the model is 
    // switched on, this function will need minor modifications.

    // Gas constant [J / mol K] 
    const scalar R = 8.314;

    //  Strength coeff. [Mpa^-n * s^-1]
    scalar A(GREAT);
    //  Activation energy for creep deformation [J/mol]
    scalar Q(GREAT);
    //  Stress exponent [-]
    scalar n(GREAT);

    // Useful quantities for linear interpolations:
    const scalar lnAalpha = 9.07532216; // == ln(8737)
    const scalar lnAbeta  = 2.06686276; // == ln(7.9)
    const scalar Qalpha   = 3.21e5 + 24.69 * (T - 923.15);
    const scalar Qbeta    = 1.41919e5;
    const scalar nAlpha   = 5.89;
    const scalar nBeta    = 3.78;

    // Initialize creep thermal rate in 1/s
    scalar creepThRate(creepThRateOld/3600);

    // Tolerance and error
    const scalar tol(1);
    scalar       err(1);

    // The following do-while loop is necessary to determine A,Q and n in case
    //  of mixed-phase Zirconium.
    do
    {
        // Initialize trial creepThRate in 1/s
        scalar creepThRateTrial(creepThRate);
        
        // Alpha Phase Zr:
        if( betaFrac == 0 )
        {
            A = 8737;
            Q = 3.21e5 + 24.69 * (T - 923.15);
            n = 5.89;
        }
        // Beta Phase Zr:
        else if ( betaFrac > 0.9999 )
        {
            A = 7.9;
            Q = 1.41919e5;
            n = 3.78;
        }
        // Mixed phase, creep strain rate <= 3e-3 (s^-1)
        else if( creepThRateTrial <= 3e-3 )
        {
            // Useful quantities for the interpolation:
            const scalar Amix   = 0.24;
            const scalar Qmix   = 1.02366e5;
            const scalar nMix   = 2.33;

            // Depending if beta is < or >= 0.5, linear interpolation
            //  of ln(A), Q and n is performed between pure alpha and
            //  50%/50% mixture or viceversa.
            if (betaFrac < 0.5)
            {
                A = exp( lnAalpha + (betaFrac-0)/(0.5-0) * (log(Amix)-lnAalpha) );
                Q = Qalpha + (betaFrac-0)/(0.5-0) * (Qmix-Qalpha);
                n = nAlpha + (betaFrac-0)/(0.5-0) * (nMix-nAlpha);
            }
            else
            {
                A = exp( log(Amix) + (betaFrac-0.5)/(1.0-0.5) * (lnAbeta-log(Amix)) );
                Q = Qmix + (betaFrac-0.5)/(1.0-0.5) * (Qbeta-Qmix);
                n = nMix + (betaFrac-0.5)/(1.0-0.5) * (nBeta-nMix);
            }
        }
        // Mixed phase, creep strain rate > 3e-3 (s^-1)
        else if( creepThRateTrial > 3e-3 )
        {
            // In this case instead it is assumed that ln(AA), nn and QQ vary linearly 
            // between the values for pure alpha and pure beta phase.
            A = exp( lnAalpha + betaFrac * (lnAbeta-lnAalpha) );
            Q = Qalpha + betaFrac * (Qbeta-Qalpha);
            n = nAlpha + betaFrac * (nBeta-nAlpha);
        }


        // Norton power equation for Cladding large strain rate:
        creepThRate = 
        A * exp(-Q / (R*T)) * pow(sigmaEff/1e6,n) ;

        // - Update error
        err = mag(creepThRate - creepThRateTrial);

    }while( err > tol );

    // Info << "LOCA creep rate : " << creepThRate << endl;

    // Return creep rate in 1/hr
    return creepThRate* 3600;
}


Foam::scalar Foam::LimbackCreepModelLOCA::calcLOCAThermalCreepRatePrime
(
    const scalar& T,
    const scalar& betaFrac,
    const scalar& sigmaEff,
    const scalar& creepRate
)
{
    // TODO (EB) : this model makes use of the betaFraction field. Currently 
    // OFFBEAT creates the field (and initializes it to zero) even if the model
    // is explicitly set to "none" (coeherently to what swellingModel does). 
    // Now this is not a problem, as in this case the correlations will be 
    // evaluated as the material was pure alpha. In the future, if it is
    // decided that the field betaFraction is created only when the model is 
    // switched on, this function will need minor modifications.

    //  Stress exponent [-]
    scalar n(GREAT);

    const scalar nAlpha   = 5.89;
    const scalar nBeta    = 3.78;
        
    // Alpha Phase Zr:
    if( betaFrac == 0 )
    {
        n = 5.89;
    }
    // Beta Phase Zr:
    else if ( betaFrac > 0.9999 )
    {
        n = 3.78;
    }
    // Mixed phase, creep strain rate <= 3e-3 (s^-1)
    else if( creepRate/3600 <= 3e-3 )
    {
        // Useful quantities for the interpolation:
        const scalar nMix   = 2.33;

        // Depending if beta is < or >= 0.5, linear interpolation
        //  of ln(A), Q and n is performed between pure alpha and
        //  50%/50% mixture or viceversa.
        if (betaFrac < 0.5)
        {
            n = nAlpha + (betaFrac-0)/(0.5-0) * (nMix-nAlpha);
        }
        else
        {
            n = nMix + (betaFrac-0.5)/(1.0-0.5) * (nBeta-nMix);
        }
    }
    // Mixed phase, creep strain rate > 3e-3 (s^-1)
    else if( creepRate/3600  > 3e-3 )
    {
        n = nAlpha + betaFrac * (nBeta-nAlpha);
    }

    // Norton power equation for Cladding large strain rate:
    return n*creepRate/(sigmaEff/1e6) ;
}

// ************************************************************************* //