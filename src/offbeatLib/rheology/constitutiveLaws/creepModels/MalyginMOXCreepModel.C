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

#include "MalyginMOXCreepModel.H"
#include "addToRunTimeSelectionTable.H"
#include "globalFieldLists.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * *//

namespace Foam
{
    defineTypeNameAndDebug(MalyginMOXCreepModel, 0);
    addToRunTimeSelectionTable
    (
        creepModel, 
        MalyginMOXCreepModel, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * ** * * * * //

Foam::MalyginMOXCreepModel::MalyginMOXCreepModel
(
    const fvMesh& mesh,
    const dictionary& lawDict
)
:
    creepModel(mesh, lawDict),
    epsilonCreepIrr_(createOrLookup<symmTensor>(mesh, "epsilonCreepIrr")),
    epsilonCreepTh_(createOrLookup<symmTensor>(mesh, "epsilonCreepTh")),
    epsilonCreepIrrEq_(createOrLookup<scalar>(mesh, "epsilonCreepIrrEq")),
    DepsilonCreepIrrEq_(createOrLookup<scalar>(mesh, "DepsilonCreepIrrEq")),
    epsilonCreepThEq_(createOrLookup<scalar>(mesh, "epsilonCreepThEq")),
    DepsilonCreepThEq_(createOrLookup<scalar>(mesh, "DepsilonCreepThEq")),
    heatSourceName_(lawDict.lookupOrDefault<word>("heatSourceName", "Q" )),
    grainRadiusName_(lawDict.lookupOrDefault<word>("grainRadiusName", "grainRadius" )),
    oxygenMetalRatioName_(lawDict.lookupOrDefault<word>("oxygenMetalRatioName", "oxygenMetalRatio" )),
    youngModulusName_(lawDict.lookupOrDefault<word>("youngModulusName", "E" )),
    porosityName_(lawDict.lookupOrDefault<word>("porosityName", "porosity" )),
    PuConc_(readScalar(lawDict.lookup("PuConcentration"))),
    Efission_(3.204e-11),
    KB_(8.6173303e-5),
    DH_(2.5),
    Q_(0.9),
    Q1_(1.6),
    Q2_(6.9),
    Q3_(3.8),
    A1_(3.8e-14),
    A2_(2.3e-11),
    A3_(3.1e10),
    B1_(3.67e7),
    irradiationCreep_(lawDict.lookupOrDefault<bool>("irradiationCreep", true )),
    timeIndex_(-1),
    relax_(lawDict.lookupOrDefault<scalar>("relax", 1.0 )),
    NewtonRaphsonMethod_(lawDict.lookupOrDefault<bool>("NewtonRaphsonMethod", false)),
    NewtonRaphsonTolerance_(NewtonRaphsonMethod_ ? 
        lawDict.lookupOrDefault<scalar>("NewtonRaphsonTolerance", 1e-6) : 1.0),
    creepRateLimiter_(lawDict.lookupOrDefault<bool>("creepRateLimiter", false )),
    creepRateLimit_(lawDict.lookupOrDefault<scalar>("creepRateLimit", 4.1e-6 ))
{
    epsilonCreepIrr_.oldTime();
    epsilonCreepTh_.oldTime();
    epsilonCreepIrrEq_.oldTime();
    DepsilonCreepIrrEq_.oldTime();
    epsilonCreepThEq_.oldTime();
    DepsilonCreepThEq_.oldTime();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::MalyginMOXCreepModel::~MalyginMOXCreepModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::MalyginMOXCreepModel::correctCreep
(
    volSymmTensorField& epsilonEl,
    const labelList& addr
)
{
    const volScalarField& mu(mesh_.lookupObject<volScalarField>("mu"));
    const volScalarField& lambda(mesh_.lookupObject<volScalarField>("lambda"));
    const volScalarField& Temp(mesh_.lookupObject<volScalarField>("T"));

    const volScalarField& heatSource(mesh_.lookupObject<volScalarField>(heatSourceName_));
    const volScalarField& grainRadius(mesh_.lookupObject<volScalarField>(grainRadiusName_));
    const volScalarField& oxygenMetalRatio = mesh_.lookupObject<volScalarField>(oxygenMetalRatioName_);    
    const volScalarField& youngModulus = mesh_.lookupObject<volScalarField>(youngModulusName_);
    const volScalarField& porosity = mesh_.lookupObject<volScalarField>(porosityName_);

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

        const scalarField& epsilonCreepEqOld_ = 
        epsilonCreepEq_.oldTime().internalField();

        const symmTensorField& epsilonCreepIrrOldI = 
        epsilonCreepIrr_.oldTime().internalField();
        const symmTensorField& epsilonCreepThOldI = 
        epsilonCreepTh_.oldTime().internalField();
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

                            scalarField& creepThermalRateP = 
                            creepThermalRate.boundaryFieldRef()[patchID];

                            DepsilonCreepEqP[faceID] *= 0;
                            DepsilonCreepIrrEqP[faceID] *= 0;
                            DepsilonCreepThP[faceID] *= 0;

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

                // --- Useful quantities

                // Reference to temperature field [K]
                const scalar T = Temp.internalField()[cellI];
                
                maxTemp = max(T, maxTemp);

                // Fission rate in 1 / cm3s
                scalar F = heatSource[cellI]/Efission_/1e6;

                // Young modulus in Mpa
                const scalar YM = youngModulus[cellI]/1e6;

                // porosity 
                const scalar poro = porosity[cellI];
                
                // Grain size in um
                scalar G = grainRadius[cellI]*2*1e6;

                // OM ratio
                const scalar OM = oxygenMetalRatio[cellI];

                // Departure from stoichiometry
                const scalar x = OM -2;

                // sigmaEff in MPa limited at 150 MPa 
                scalar sEffMPa = min(sigmaEff/1e6, 150);

                //--- Calculate irradiation induced creep
                if (irradiationCreep_)
                {
                    // Compute lambda
                    scalar L = 1e11 * exp(-min(4.5/(KB_*T), 70));

                    // Compute Kvs: formation rate of mobile radiation defects 
                    scalar Kvs = (T<1173)
                    ? F/101
                    : F/(1+(1e-18 * F)/(1e-20 * F + L));

                    // Compute  irradiation creep rate
                    scalar comp1 = 
                    A1_ * sEffMPa * Kvs / YM;

                    scalar comp2 = 
                    A2_ * sEffMPa * pow(Kvs, 0.5) * exp(-min(Q_/(KB_*T), 70));

                    scalar comp3 = 
                    1 + 0.31*pow(poro*100, 1.8);

                    creepIrrRate = (comp1 + comp2) * comp3;

                    if (NewtonRaphsonMethod_)
                    {
                        // Compute derivative wrt stress
                        creepIrrRatePrime = creepIrrRate/sEffMPa;
                    }

                    // Update irradiation creep equivalent
                    DepsilonCreepIrrEq_.ref()[cellI] =  
                    (1 - relax_)*DepsilonCreepIrrEq_.ref()[cellI] 
                    + relax_*creepIrrRate*deltaT;

                    epsilonCreepIrrEq_.ref()[cellI] = 
                    epsilonCreepIrrEqOld_[cellI] + DepsilonCreepIrrEq_.ref()[cellI];
                }

                //--- Calculate thermal creep (eq5 of reference vol 1)
                
                scalar R = pow(sqrt(x*x + 82*exp(-min(DH_/(KB_*T),70))) - x, 2.0);

                scalar A = A3_*exp(19.0*poro + 3.6*PuConc_);
                scalar B = B1_*exp(10.5*poro + 3.6*PuConc_);
                
                scalar factor1 = A * sEffMPa / pow(G, 2.0)
                               + B * pow(sEffMPa, 4.5) * exp(-min(Q1_/(KB_*T), 70));

                scalar factor2 = 165/R * exp(-min(Q2_/(KB_*T), 70))
                               + R/(1+450*x*x)* exp(-min(Q3_/(KB_*T), 70));

                scalar creepThRate = ( T >=1273.15 ) 
                ? factor1 * factor2
                : 0.0;

                if (NewtonRaphsonMethod_)
                {
                    // Compute derivative of thermal creep wrt sigmaEff
                    creepThRatePrime = (A / pow(G, 2.0) + 4.5 * B / sEffMPa * pow(sEffMPa, 4.5) * exp(-min(Q1_/(KB_*T), 70))) * factor2;
                }

                creepThermalRate[cellI] = creepThRate;

                // Update thermal creep equivalent
                DepsilonCreepThEq_.ref()[cellI] = 
                (1 - relax_)*DepsilonCreepThEq_.ref()[cellI] + 
                relax_*creepThRate*deltaT;

                epsilonCreepThEq_.ref()[cellI] = 
                epsilonCreepThEqOld_[cellI] + DepsilonCreepThEq_.ref()[cellI];

                // --- Total creep rate

                scalar fTrial
                (
                    creepIrrRate 
                  + creepThRate 
                );

                if ( creepRateLimiter_ )
                {
                    // Limit Creep rate 
                    fTrial = min(fTrial, creepRateLimit_);
                }

                if (NewtonRaphsonMethod_)
                {
                    const scalar fTrialPrime
                    (
                        creepIrrRatePrime
                      + creepThRatePrime 
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
                            const symmTensorField& epsilonElP = 
                            epsilonEl.boundaryField()[patchID];

                            const scalarField& grainRadiusP =
                            grainRadius.boundaryField()[patchID];
                            const scalarField& heatSourceP =
                            heatSource.boundaryField()[patchID];
                            const scalarField& youngModulusP = 
                            youngModulus.boundaryField()[patchID];
                            const scalarField& porosityP = 
                            porosity.boundaryField()[patchID];
                            const scalarField& oxygenMetalRatioP = 
                            oxygenMetalRatio.boundaryField()[patchID];
                           


                            symmTensorField& epsilonCreepP = 
                            epsilonCreep_.boundaryFieldRef()[patchID];
                            symmTensorField& DepsilonCreepP = 
                            DepsilonCreep_.boundaryFieldRef()[patchID];
                            symmTensorField& epsilonCreepIrrP = 
                            epsilonCreepIrr_.boundaryFieldRef()[patchID];
                            symmTensorField& epsilonCreepThP = 
                            epsilonCreepTh_.boundaryFieldRef()[patchID];

                            scalarField& epsilonCreepEqP = 
                            epsilonCreepEq_.boundaryFieldRef()[patchID];
                            scalarField& epsilonCreepIrrEqP = 
                            epsilonCreepIrrEq_.boundaryFieldRef()[patchID];
                            scalarField& epsilonCreepThEqP = 
                            epsilonCreepThEq_.boundaryFieldRef()[patchID];

                            const scalarField& epsilonCreepIrrEqOldP = 
                            epsilonCreepIrrEq_.oldTime().boundaryField()[patchID];
                            const scalarField& epsilonCreepThEqOldP = 
                            epsilonCreepThEq_.oldTime().boundaryField()[patchID];
                            const scalarField& epsilonCreepEqOldP = 
                            epsilonCreepEq_.oldTime().boundaryField()[patchID];

                            const symmTensorField& epsilonCreepIrrOldP = 
                            epsilonCreepIrr_.oldTime().boundaryField()[patchID];
                            const symmTensorField& epsilonCreepThOldP = 
                            epsilonCreepTh_.oldTime().boundaryField()[patchID];
                            const symmTensorField& epsilonCreepOldP = 
                            epsilonCreep_.oldTime().boundaryField()[patchID];

                            scalarField& DepsilonCreepIrrEqP = 
                            DepsilonCreepIrrEq_.boundaryFieldRef()[patchID];
                            scalarField& DepsilonCreepThEqP = 
                            DepsilonCreepThEq_.boundaryFieldRef()[patchID];
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

                            // --- Useful quantities

                            // Reference to temperature field [K]
                            const scalar T = Tp[faceID];
                            
                            // Fission rate in 1 / cm3s
                            scalar F = heatSourceP[faceID]/Efission_/1e6;

                            // Young modulus in Mpa
                            const scalar YM = youngModulusP[faceID]/1e6;

                            // porosity 
                            const scalar poro = porosityP[faceID];
                            
                            // Grain size in um
                            scalar G = grainRadiusP[faceID]*2*1e6;

                            // OM ratio
                            const scalar OM = oxygenMetalRatioP[faceID];

                            // Departure from stoichiometry
                            const scalar x = OM -2;

                            // sigmaEff in MPa limited at 150 MPa 
                            scalar sEffMPa = min(sigmaEff/1e6, 150);

                            //--- Calculate irradiation induced creep

                            if (irradiationCreep_)
                            {
                                // Compute lambda
                                scalar L = 1e11 * exp(-min(4.5/(KB_*T), 70));

                                // Compute Kvs: formation rate of mobile radiation defects 
                                scalar Kvs = (T<1173)
                                ? F/101
                                : F/(1+(1e-18 * F)/(1e-20 * F + L));

                                // Compute  irradiation creep rate
                                scalar comp1 = 
                                A1_ * sEffMPa * Kvs / YM;

                                scalar comp2 = 
                                A2_ * sEffMPa * pow(Kvs, 0.5) * exp(-min(Q_/(KB_*T), 70));

                                scalar comp3 = 
                                1 + 0.31*pow(poro*100, 1.8);

                                creepIrrRate = (comp1 + comp2) * comp3;

                                if (NewtonRaphsonMethod_)
                                {
                                    // Compute derivative wrt stress
                                    creepIrrRatePrime = creepIrrRate/sEffMPa;
                                }

                                DepsilonCreepIrrEqP[faceID] =  
                                (1 - relax_)*DepsilonCreepIrrEqP[faceID] 
                                + relax_*creepIrrRate*deltaT;

                                epsilonCreepIrrEqP[faceID] = 
                                DepsilonCreepIrrEqP[faceID] + epsilonCreepIrrEqOldP[faceID];
                            }

                            //--- Calculate thermal creep (eq5 of reference vol 1)
            
                            scalar R = pow(sqrt(x*x + 82*exp(-min(DH_/(KB_*T),70))) - x, 2.0);

                            scalar A = A3_*exp(19.0*poro + 3.6*PuConc_);
                            scalar B = B1_*exp(10.5*poro + 3.6*PuConc_);
                            
                            scalar factor1 = A * sEffMPa / pow(G, 2.0)
                                        + B * pow(sEffMPa, 4.5) * exp(-min(Q1_/(KB_*T), 70));

                            scalar factor2 = 165/R * exp(-min(Q2_/(KB_*T), 70))
                                        + R/(1+450*x*x)* exp(-min(Q3_/(KB_*T), 70));

                            scalar creepThRate = ( T >=1273.15 ) 
                            ? factor1 * factor2
                            : 0.0;

                            if (NewtonRaphsonMethod_)
                            {
                                // Compute derivative of thermal creep wrt sigmaEff
                                creepThRatePrime = (A / pow(G, 2.0) + 4.5 * B / sEffMPa * pow(sEffMPa, 4.5) * exp(-min(Q1_/(KB_*T), 70))) * factor2;
                            }

                            creepThermalRateP[faceID] = creepThRate;

                            // Update thermal creep equivalent
                            DepsilonCreepThEqP[faceID] = 
                            (1 - relax_)*DepsilonCreepThEqP[faceID] + 
                            relax_*creepThRate*deltaT;               

                            epsilonCreepThEqP[faceID] = 
                            DepsilonCreepThEqP[faceID] + epsilonCreepThEqOldP[faceID] ;

                            /******* Total creep rate ********/

                            scalar fTrial
                            (
                                creepIrrRate 
                              + creepThRate
                            );

                            if ( creepRateLimiter_ )
                            {
                                // Limit Creep rate 
                                fTrial = min(fTrial, creepRateLimit_);
                            }

                            maxRate = max(fTrial, maxRate);

                            if (NewtonRaphsonMethod_)
                            {
                                const scalar fTrialPrime
                                (
                                    creepIrrRatePrime
                                  + creepThRatePrime 
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

// ************************************************************************* //