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

#include "MatproCreepModel.H"
#include "addToRunTimeSelectionTable.H"
#include "wedgeFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(MatproCreepModel, 0);
    addToRunTimeSelectionTable
    (
        creepModel, 
        MatproCreepModel, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MatproCreepModel::MatproCreepModel
(
    const fvMesh& mesh,
    const dictionary& lawDict
)
:
    creepModel(mesh, lawDict),
    heatSourceName_(lawDict.lookupOrDefault<word>("heatSourceName", "Q" )),
    densityName_(lawDict.lookupOrDefault<word>("densityName", "rho" )),
    grainRadiusName_(lawDict.lookupOrDefault<word>("grainRadiusName", "grainRadius" )),
    oxygenMetalRatioName_(lawDict.lookupOrDefault<word>("oxygenMetalRatioName", "oxygenMetalRatio" )),
    SakaiCorrection_(lawDict.lookupOrDefault<Switch>("SakaiCorrection", true)),
    A1(0.3919),
    A2(1.3100e-19),
    A3(-87.7),
    A4(2.0391e-25),
    A6(-90.5),
    A7(SakaiCorrection_ ? 7.78e-37 : 3.7226e-35),
    R(8.3143),
    timeIndex_(-1),
    relax_(lawDict.lookupOrDefault<scalar>("relax", 1.0 )),
    NewtonRaphsonMethod_(lawDict.lookupOrDefault<bool>("NewtonRaphsonMethod", false)),
    NewtonRaphsonTolerance_(NewtonRaphsonMethod_ ? 
        lawDict.lookupOrDefault<scalar>("NewtonRaphsonTolerance", 1e-6) : 1.0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::MatproCreepModel::~MatproCreepModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::MatproCreepModel::correctCreep
(
    volSymmTensorField& epsilonEl,
    const labelList& addr
)
{ 
    // epsilonCreep_.storeOldTimes();
    // epsilonCreepEq_.storeOldTimes();

    const volScalarField& mu(mesh_.lookupObject<volScalarField>("mu"));
    const volScalarField& lambda(mesh_.lookupObject<volScalarField>("lambda"));
    const volScalarField& Temp(mesh_.lookupObject<volScalarField>("T"));

    scalar tol(1.0);
    //- Irradiation time step [s]
    scalar deltaT(0);
    
    if(timeIndex_ > 0)
    {
        const volScalarField& heatSource = mesh_.lookupObject<volScalarField>(heatSourceName_);
        const volScalarField& grainRadius = mesh_.lookupObject<volScalarField>(grainRadiusName_);
        const volScalarField& rho = mesh_.lookupObject<volScalarField>(densityName_);    
        const volScalarField& oxygenMetalRatio = mesh_.lookupObject<volScalarField>(oxygenMetalRatioName_);    

        deltaT = mesh_.time().deltaTValue();

        const scalarField& epsilonCreepEqOld_ = epsilonCreepEq_.oldTime().internalField();
        const symmTensorField& epsilonCreepOldI = epsilonCreep_.oldTime().internalField();
    
        //- Adjust elastic strain tensor by subtracting old creep values    
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

        // if(mesh_.time().timeIndex() > timeIndex_)
        // {
        //     forAll(addr, addrI)
        //     {
        //         const label cellI = addr[addrI];
        //         DepsilonCreepEq_.ref()[cellI] *= 0;

        //         const cell& c = mesh_.cells()[cellI];  

        //         forAll(c, faceI)
        //         {
        //             const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);

        //             {
        //                 if (patchID > -1 and epsilonEl.boundaryField()[patchID].size())
        //                 {
        //                     const label faceID = 
        //                     mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);

        //                     scalarField& DepsilonCreepEqP = DepsilonCreepEq_.boundaryFieldRef()[patchID];

        //                     DepsilonCreepEqP[faceID] *= 0;
        //                 }
        //             }
        //         }
        //     }
        // }

        do
        {
            scalar averageCreepEq(0);
            scalar averageRate(0);
            scalar maxRate(0);
            scalar maxDeltaCreep(0);
            scalar maxFissionRate(0);
            scalar maxTemp(0);
            scalar maxStress(0);
            scalar volTot(0);

            
            scalar creepIncrementPrime1(0);
            scalar creepIncrementPrime2(0);
            scalar creepIncrementPrime3(0);
            scalar fTrialDev(0);
            scalar deltaCreep(0);

            forAll(addr, addrI)
            {
                const label cellI = addr[addrI];

                const symmTensor sigmaI = 
                2*mu[cellI]*epsilonEl[cellI] + lambda[cellI]*tr(epsilonEl[cellI])*I;

                //- Calculate the von Mises stress [Pa]
                const scalar sigmaEff = 
                max
                (
                    min
                    (
                        sqrt((3.0/2.0)*magSqr(dev(sigmaI))) 
                        - 3*mu[cellI]*DepsilonCreepEq_.ref()[cellI],
                        1e10
                    ), 
                    1e5
                );

                const scalar DepsilonCreepEqPrev = DepsilonCreepEq_.ref()[cellI];

                //- Reference to temperature field [K]
                const scalar T = Temp.internalField()[cellI];

                /****** Calculate activation energies and other parameters  *********/  
                const scalar OM = oxygenMetalRatio[cellI];

                scalar f = 1/(exp(-20/Foam::log(max(OM - 2, 1e-15)) - 8) + 1);

                scalar Q1 = 74829*f + 301762;
                scalar Q2 = 83143*f + 469191;
                scalar Q3 = 21759;

                //- Fission rate
                scalar F = heatSource[cellI]/312e-13;
                //- Grain size in um
                scalar G = grainRadius[cellI]*2*1e6;
                //- Density Fraction
                scalar D = rho[cellI]/10970*100;

                // scalar sigmaTransition = 1.6547e7/pow(G, 0.5714);
                // scalar sigmaTransition = GREAT;

                /****** Calculate first thermal creep term *********/  

                scalar creepIncrement1 = 
                ( (A1 + A2*F)/((A3 + D)*G*G) )
                *sigmaEff*exp(-Q1/R/T);

                if (NewtonRaphsonMethod_)
                {
                    creepIncrementPrime1 = 
                    ( (A1 + A2*F)/((A3 + D)*G*G) )*exp(-Q1/R/T);
                }

                /****** Calculate second thermal creep term *********/  

                scalar creepIncrement2 = ( A4/(A6 + D) )*pow(sigmaEff, 4.5)*exp(-Q2/R/T);
                
                if (NewtonRaphsonMethod_)
                {
                    creepIncrementPrime2 = 4.5*( A4/(A6 + D) )*pow(sigmaEff, 3.5)*exp(-Q2/R/T);
                }

                // if(sigmaEff < sigmaTransition)
                // {
                //     creepIncrement2 *= 0;
                //     creepIncrementPrime2 *= 0;
                // }

                /****** Calculate third creep term  *********/  

                scalar creepIncrement3(0.0);
                if(SakaiCorrection_)
                {
                    creepIncrement3 = A7*F*sigmaEff;
                    if (NewtonRaphsonMethod_)
                    {
                        creepIncrementPrime3 = A7*F;
                    }
                }
                else
                {
                    creepIncrement3 = A7*F*sigmaEff*exp(-Q3/R/T);
                    if (NewtonRaphsonMethod_)
                    {
                        creepIncrementPrime3 = A7*F*exp(-Q3/R/T);
                    }
                }

                /****** Calculate new creep increment  *********/  

                // Info << "T: " << T << ", sigmaEff: " << sigmaEff << ", f: " << f << ", F: " << F << ", G: " << G << ", D: " << D << ", deltaT: " << deltaT << endl;
                // Info << "sigmaEff: " << sigmaEff<< ", creepIncrement1: " << creepIncrement1 << ", creepIncrement2: " << creepIncrement2 << ", creepIncrement3: " << creepIncrement3 << endl << endl;

                //- Creep and Creep derivative
                scalar fTrial = min
                (
                    (creepIncrement1 + creepIncrement2 + creepIncrement3),
                    GREAT
                    // 2.78e-9
                );
                
                if (NewtonRaphsonMethod_)
                {
                    fTrialDev = (creepIncrementPrime1 + creepIncrementPrime2 + creepIncrementPrime3);
                    deltaCreep = (DepsilonCreepEq_.ref()[cellI] - (DepsilonCreepEq_.ref()[cellI] - deltaT*fTrial)/(1 + 3*mu[cellI]*fTrialDev*deltaT));
                }
                else
                {
                    deltaCreep = fTrial*deltaT;
                }

                averageRate += fTrial*mesh_.V()[cellI];
                averageCreepEq += fTrial*deltaT*mesh_.V()[cellI];
                volTot += mesh_.V()[cellI];
                maxRate = max(fTrial, maxRate);
                maxDeltaCreep = max(deltaCreep, maxDeltaCreep);
                maxFissionRate = max(F, maxFissionRate);
                maxTemp = max(T, maxTemp);
                maxStress = max(sigmaEff, maxStress);

                /****** Update creep  *********/  
                DepsilonCreepEq_.ref()[cellI] = (1 - relax_)*DepsilonCreepEq_.ref()[cellI] + relax_*deltaCreep;   
                
                epsilonCreepEq_.ref()[cellI] = epsilonCreepEqOld_[cellI] + DepsilonCreepEq_.ref()[cellI];

                DepsilonCreep_.ref()[cellI] = (1 - relax_)*DepsilonCreep_.ref()[cellI] + relax_*1.5*DepsilonCreepEq_.ref()[cellI]*dev(sigmaI - 2*mu[cellI]*DepsilonCreep_.ref()[cellI])/max(sigmaEff, VSMALL); 

                epsilonCreep_.ref()[cellI] = epsilonCreepOldI[cellI] + DepsilonCreep_.ref()[cellI]; 

                const scalar tolTr = min(1, mag( (DepsilonCreepEq_.ref()[cellI] - DepsilonCreepEqPrev)/ (max(DepsilonCreepEqPrev , VSMALL))));

                tol = min(tol, tolTr );                

                const cell& c = mesh_.cells()[cellI];  

                forAll(c, faceI)
                {
                    const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);

                    {
                        if (patchID > -1 and epsilonEl.boundaryField()[patchID].size() and 
                            not(isType<wedgeFvPatchField<symmTensor>>(epsilonEl.boundaryField()[patchID])))
                        {
                            const label faceID = 
                            mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);

                            const scalarField& Tp = Temp.boundaryField()[patchID];
                            const scalarField& heatSourceP = heatSource.boundaryField()[patchID];
                            const scalarField& grainRadiusP = grainRadius.boundaryField()[patchID];

                            const scalarField& rhoP = rho.boundaryField()[patchID];
                            const scalarField& oxygenMetalRatioP = oxygenMetalRatio.boundaryField()[patchID];
                            const scalarField& muP = mu.boundaryField()[patchID];
                            const scalarField& lambdaP = lambda.boundaryField()[patchID];
                            const symmTensorField& epsilonElP = epsilonEl.boundaryField()[patchID];

                            symmTensorField& epsilonCreepP = epsilonCreep_.boundaryFieldRef()[patchID];
                            symmTensorField& DepsilonCreepP = DepsilonCreep_.boundaryFieldRef()[patchID];
                            const symmTensorField& epsilonCreepOldP = epsilonCreep_.oldTime().boundaryField()[patchID];
                            
                            scalarField& epsilonCreepEqP = epsilonCreepEq_.boundaryFieldRef()[patchID];
                            scalarField& DepsilonCreepEqP = DepsilonCreepEq_.boundaryFieldRef()[patchID];
                            const scalarField& epsilonCreepEqOldP = epsilonCreepEq_.oldTime().boundaryField()[patchID];  

                            //- Calculate the von Mises stress [Pa]
                            const symmTensor sigmaPI = 
                            2*muP[faceID]*epsilonElP[faceID] + lambdaP[faceID]*tr(epsilonElP[faceID])*I;

                            const scalar sigmaEff = max(min(sqrt((3.0/2.0)*magSqr(dev(sigmaPI))) - 3*muP[faceID]*DepsilonCreepEqP[faceID], 1e10), 1e5);

                            const scalar DepsilonCreepEqPPrev = DepsilonCreepEqP[faceID];

                            //- Reference to temperature field [K]
                            const scalar T = Tp[faceID];

                             /****** Calculate activation energies and other parameters  *********/  
                            const scalar OM = oxygenMetalRatioP[faceID];

                            scalar f = 1/(exp(-20/Foam::log(max(OM - 2, 1e-15)) - 8) + 1);

                            scalar Q1 = 74829*f + 301762;
                            scalar Q2 = 83143*f + 469191;
                            scalar Q3 = 21759;

                            //- Fission rate
                            scalar F = heatSourceP[faceID]/312e-13;
                            //- Grain size in um
                            scalar G = grainRadiusP[faceID]*2*1e6;
                            //- Density Fraction
                            scalar D = rhoP[faceID]/10970*100;

                            // scalar sigmaTransition = 1.6547e7/pow(G, 0.5714);
                            // scalar sigmaTransition = GREAT;

                            /****** Calculate first thermal creep term *********/   

                            scalar creepIncrement1 = 
                            ( (A1 + A2*F)/((A3 + D)*G*G) )
                            *sigmaEff*exp(-Q1/R/T);

                            if (NewtonRaphsonMethod_)
                            {
                                creepIncrementPrime1 = ( (A1 + A2*F)/((A3 + D)*G*G) )*exp(-Q1/R/T);
                            }

                            /****** Calculate second thermal creep term *********/  

                            scalar creepIncrement2 = ( A4/(A6 + D) )*pow(sigmaEff, 4.5)*exp(-Q2/R/T);
                            if (NewtonRaphsonMethod_)
                            {
                                creepIncrementPrime2 = 4.5*( A4/(A6 + D) )*pow(sigmaEff, 3.5)*exp(-Q2/R/T);
                            }
                            
                            // if(sigmaEff < sigmaTransition)
                            // {
                            //     creepIncrement2 *= 0;
                            //     creepIncrementPrime2 *= 0;
                            // }

                            /****** Calculate third creep term  *********/  

                            scalar creepIncrement3(0.0);
                            if(SakaiCorrection_)
                            {
                                creepIncrement3 = A7*F*sigmaEff;
                                if (NewtonRaphsonMethod_)
                                {
                                    creepIncrementPrime3 = A7*F;
                                }
                            }
                            else
                            {
                                creepIncrement3 = A7*F*sigmaEff*exp(-Q3/R/T);
                                if (NewtonRaphsonMethod_)
                                {
                                    creepIncrementPrime3 = A7*F*exp(-Q3/R/T);
                                }
                            }

                            /****** Calculate new creep increment  *********/  
                            
                            // Info << "FACE: sigmaEff: " << sigmaEff<< ", creepIncrement1: " << creepIncrement1 << ", creepIncrement2: " << creepIncrement2 << ", creepIncrement3: " << creepIncrement3 << endl << endl;

                            //- Creep and Creep derivative
                            scalar fTrial = min
                            (
                                (creepIncrement1 + creepIncrement2 + creepIncrement3),
                                GREAT
                                // 2.78e-9
                            );

                            if (NewtonRaphsonMethod_)
                            {
                                fTrialDev = (creepIncrementPrime1 + creepIncrementPrime2 + creepIncrementPrime3);
                                deltaCreep = (DepsilonCreepEqP[faceID] - (DepsilonCreepEqP[faceID] - deltaT*fTrial)/(1 + 3*mu[faceID]*fTrialDev*deltaT));
                            }
                            else
                            {
                                deltaCreep = fTrial*deltaT;
                            }

                            //- New increment with Newton Ralphson method
                            // maxRate = max(fTrial, maxRate);
                            // maxFissionRate = max(F, maxFissionRate);
                            // maxTemp = max(T, maxTemp);
                            // maxStress = max(sigmaEff, maxStress);


                        /****** Update creep  *********/  
                            DepsilonCreepEqP[faceID] = (1-relax_)*DepsilonCreepEqP[faceID] + relax_*deltaCreep;  
                                
                            epsilonCreepEqP[faceID] = DepsilonCreepEqP[faceID] + epsilonCreepEqOldP[faceID] ;

                            DepsilonCreepP[faceID] = (1 - relax_)*DepsilonCreepP[faceID] + relax_*1.5*DepsilonCreepEqP[faceID]*dev(sigmaPI - 2*muP[faceID]*DepsilonCreepP[faceID])/max(sigmaEff, VSMALL) ;
                            
                            epsilonCreepP[faceID] = epsilonCreepOldP[faceID] + DepsilonCreepP[faceID];
                
                            const scalar tolTr = min(1, mag( (DepsilonCreepEqP[faceID] - DepsilonCreepEqPPrev)/ (max(DepsilonCreepEqPPrev , VSMALL))));

                            tol = min(tol, tolTr );
                        }
                    }
                }
            }

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
    
    timeIndex_ = mesh_.time().timeIndex();  

}

// ************************************************************************* //


