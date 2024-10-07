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

#include "TobbeDINCreepModel.H"
#include "addToRunTimeSelectionTable.H"
#include "globalFieldLists.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(TobbeDINCreepModel, 0);
    addToRunTimeSelectionTable
    (
        creepModel, 
        TobbeDINCreepModel, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::TobbeDINCreepModel::TobbeDINCreepModel
(
    const fvMesh& mesh,
    const dictionary& lawDict
)
:
    creepModel(mesh, lawDict),
    fastFluxName_(lawDict.lookupOrDefault<word>("fastFluxName", "fastFlux" )),
    fastFlux_(0),
    epsilonCreepIrr_(createOrLookup<symmTensor>(mesh, "epsilonCreepIrr")),
    epsilonCreepTh_(createOrLookup<symmTensor>(mesh, "epsilonCreepTh")),
    epsilonCreepIrrEq_(createOrLookup<scalar>(mesh, "epsilonCreepIrrEq")),
    DepsilonCreepIrrEq_(createOrLookup<scalar>(mesh, "DepsilonCreepIrrEq")),
    epsilonCreepThEq_(createOrLookup<scalar>(mesh, "epsilonCreepThEq")),
    DepsilonCreepThEq_(createOrLookup<scalar>(mesh, "DepsilonCreepThEq")),
    timeIndex_(0),
    relax_(lawDict.lookupOrDefault<scalar>("relax", 1.0 )),
    par1(3.2e-24),
    // par2(7.49e15),
    par2(7.52e15),
    par3(9.1e4),
    par4(1.986),
    // par5(57.46)
    par5(49.67)

{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::TobbeDINCreepModel::~TobbeDINCreepModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::TobbeDINCreepModel::correctCreep
(
    volSymmTensorField& epsilonEl,
    const labelList& addr
)
{    
    if(fastFlux_ == nullptr)
    {        
        fastFlux_    = &mesh_.lookupObject<volScalarField>(fastFluxName_);
    }   

    const volScalarField& mu(mesh_.lookupObject<volScalarField>("mu"));
    const volScalarField& lambda(mesh_.lookupObject<volScalarField>("lambda"));
    const volScalarField& Temp(mesh_.lookupObject<volScalarField>("T"));

    scalar tol(1.0);

    //- Irradiation time step [hr]
    scalar deltaT(0);
    
    if(mesh_.time().timeIndex() > 0)
    {
        deltaT = mesh_.time().deltaTValue()/3600.0;

        const scalarField& epsilonCreepIrrEqOld_ = epsilonCreepIrrEq_.oldTime().internalField();
        const scalarField& epsilonCreepThEqOld_ = epsilonCreepThEq_.oldTime().internalField();

        const scalarField& fastFluxI_ = fastFlux_->internalField();

        const scalarField& epsilonCreepEqOld_ = epsilonCreepEq_.oldTime().internalField();

        const symmTensorField& epsilonCreepIrrOldI = epsilonCreepIrr_.oldTime().internalField();
        const symmTensorField& epsilonCreepThOldI = epsilonCreepTh_.oldTime().internalField();
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

        if(mesh_.time().timeIndex() > timeIndex_)
        { 
            forAll(addr, addrI)
            {
                const label cellI = addr[addrI];
                DepsilonCreepEq_.ref()[cellI] *= 0;
                DepsilonCreepIrrEq_.ref()[cellI] *= 0;
                DepsilonCreepThEq_.ref()[cellI] *= 0;

                const cell& c = mesh_.cells()[cellI];  

                forAll(c, faceI)
                {
                    const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);

                    {
                        if (patchID > -1 and epsilonEl.boundaryField()[patchID].size())
                        {
                            const label faceID = 
                            mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);

                            scalarField& DepsilonCreepEqP = DepsilonCreepEq_.boundaryFieldRef()[patchID];
                            scalarField& DepsilonCreepIrrEqP = DepsilonCreepIrrEq_.boundaryFieldRef()[patchID];
                            scalarField& DepsilonCreepThP = DepsilonCreepThEq_.boundaryFieldRef()[patchID];

                            DepsilonCreepEqP[faceID] *= 0;
                            DepsilonCreepIrrEqP[faceID] *= 0;
                            DepsilonCreepThP[faceID] *= 0;
                        }
                    }
                }
            }

            timeIndex_ = mesh_.time().timeIndex();
        }

        do
        {
            scalar maxRate(0);
            scalar maxFlux(0);
            scalar maxTemp(0);
            scalar maxStress(0);

            forAll(addr, addrI)
            {
                const label cellI = addr[addrI];

                const symmTensor sigmaI = 
                2*mu[cellI]*epsilonEl[cellI] + lambda[cellI]*tr(epsilonEl[cellI])*I;

                //- Calculate the von Mises stress [Pa]
                const scalar sigmaEff = 
                min
                (
                    sqrt((3.0/2.0)
                    *magSqr(dev(sigmaI))) 
                    - 3*mu[cellI]*DepsilonCreepEq_.ref()[cellI] 
                    , 3e10
                );
                
                maxStress = max(sigmaEff, maxStress);

                const scalar DepsilonCreepEqPrev = DepsilonCreepEq_.ref()[cellI];

                //- Reference to temperature field [K]
                const scalar T = Temp.internalField()[cellI];

                maxTemp = max(T, maxTemp);

                /*--- IRRADIATION CREEP RATE [1/hr] ---*/ 

                //- Fast flux 
                scalar phi = fastFluxI_[cellI]; 
                maxFlux = max(phi, maxFlux);

                //- TODO (E.B.) :
                //- This correlation requires the average energy of the fast
                //- neutron flux. Here it is hardcoded to be 1 MeV.

                //- Average
                const scalar avgE = 1;

                //- This correlation requires the total neutron flux, but as 
                //- this cladding is for fast reactors, we consider the fast
                //- flux only.

                //- Tobbe correlation for 15-15 Ti austenitic [1/hr]
                const scalar creepIrrRate = 
                1e-2 * par1 * avgE * phi * sigmaEff/1e6 ;

                DepsilonCreepIrrEq_.ref()[cellI] =  
                (1 - relax_)*DepsilonCreepIrrEq_.ref()[cellI] + relax_*creepIrrRate*deltaT;

                /*--- SECONDARY THERMAL CREEP RATE [1/hr] ---*/ 

                //- Tobbe correlation for 15-15 Ti austenitic - 1.4970m
                scalar creepThRate = 1e-2 *
                par2 * exp(-par3/(par4*T)) * Foam::sinh(par5*sigmaEff/1e6/(par4*T));

                DepsilonCreepThEq_.ref()[cellI] =   
                (1 - relax_)*DepsilonCreepThEq_.ref()[cellI] + relax_*creepThRate*deltaT;

                /*--- TOTAL CREEP [1/hr] ---*/ 
                const scalar fTrial =  creepIrrRate + creepThRate;
                
                DepsilonCreepEq_.ref()[cellI] = 
                (1 - relax_)*DepsilonCreepEq_.ref()[cellI] + relax_*fTrial*deltaT;
                
                maxRate = max(fTrial, maxRate);

                epsilonCreepIrrEq_.ref()[cellI] = 
                epsilonCreepIrrEqOld_[cellI] + DepsilonCreepIrrEq_.ref()[cellI];

                epsilonCreepIrr_.ref()[cellI] = 
                epsilonCreepIrrOldI[cellI] + 1.5*DepsilonCreepIrrEq_.ref()[cellI]*dev(sigmaI 
                    - 2*mu[cellI]*DepsilonCreep_.ref()[cellI])/max(sigmaEff, SMALL);

                epsilonCreepThEq_.ref()[cellI] = 
                epsilonCreepThEqOld_[cellI] + DepsilonCreepThEq_.ref()[cellI];

                epsilonCreepTh_.ref()[cellI] = 
                epsilonCreepThOldI[cellI] + 1.5*DepsilonCreepThEq_.ref()[cellI]*dev(sigmaI  
                    - 2*mu[cellI]*DepsilonCreep_.ref()[cellI])/max(sigmaEff, SMALL); 

                epsilonCreepEq_.ref()[cellI] = 
                epsilonCreepEqOld_[cellI] + DepsilonCreepEq_.ref()[cellI];

                DepsilonCreep_.ref()[cellI] = 
                (1 - relax_)*DepsilonCreep_.ref()[cellI] + relax_*1.5*DepsilonCreepEq_.ref()[cellI]*dev(sigmaI 
                    - 2*mu[cellI]*DepsilonCreep_.ref()[cellI])/max(sigmaEff, SMALL); 

                epsilonCreep_.ref()[cellI] = 
                epsilonCreepOldI[cellI] + DepsilonCreep_.ref()[cellI]; 


                const scalar tolTr = min(1, mag( (DepsilonCreepEq_.ref()[cellI] - DepsilonCreepEqPrev)
                    / (max(DepsilonCreepEqPrev , SMALL))));

                tol = min(tol, tolTr );

                const cell& c = mesh_.cells()[cellI];  

                forAll(c, faceI)
                {
                    const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);

                    {
                        if (patchID > -1 and epsilonEl.boundaryField()[patchID].size())
                        {
                            const label faceID = 
                            mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);

                            const scalarField& Tp = Temp.boundaryField()[patchID];
                            const scalarField& muP = mu.boundaryField()[patchID];
                            const scalarField& lambdaP = lambda.boundaryField()[patchID];
                            const scalarField& fastFluxP = fastFlux_->boundaryField()[patchID];
                            const symmTensorField& epsilonElP = epsilonEl.boundaryField()[patchID];

                            symmTensorField& epsilonCreepP = epsilonCreep_.boundaryFieldRef()[patchID];
                            symmTensorField& DepsilonCreepP = DepsilonCreep_.boundaryFieldRef()[patchID];
                            symmTensorField& epsilonCreepIrrP = epsilonCreepIrr_.boundaryFieldRef()[patchID];
                            symmTensorField& epsilonCreepThP = epsilonCreepTh_.boundaryFieldRef()[patchID];

                            scalarField& epsilonCreepEqP = epsilonCreepEq_.boundaryFieldRef()[patchID];
                            scalarField& epsilonCreepIrrEqP = epsilonCreepIrrEq_.boundaryFieldRef()[patchID];
                            scalarField& epsilonCreepThEqP = epsilonCreepThEq_.boundaryFieldRef()[patchID];

                            const scalarField& epsilonCreepIrrEqOldP = epsilonCreepIrrEq_.oldTime().boundaryField()[patchID];
                            const scalarField& epsilonCreepThEqOldP = epsilonCreepThEq_.oldTime().boundaryField()[patchID];
                            const scalarField& epsilonCreepEqOldP = epsilonCreepEq_.oldTime().boundaryField()[patchID];

                            const symmTensorField& epsilonCreepIrrOldP = epsilonCreepIrr_.oldTime().boundaryField()[patchID];
                            const symmTensorField& epsilonCreepThOldP = epsilonCreepTh_.oldTime().boundaryField()[patchID];
                            const symmTensorField& epsilonCreepOldP = epsilonCreep_.oldTime().boundaryField()[patchID];

                            scalarField& DepsilonCreepIrrEqP = DepsilonCreepIrrEq_.boundaryFieldRef()[patchID];
                            scalarField& DepsilonCreepThEqP = DepsilonCreepThEq_.boundaryFieldRef()[patchID];
                            scalarField& DepsilonCreepEqP = DepsilonCreepEq_.boundaryFieldRef()[patchID];

                            //- Calculate the von Mises stress [Pa]
                            const symmTensor sigmaPI = 
                            2*muP[faceID]*epsilonElP[faceID] + lambdaP[faceID]*tr(epsilonElP[faceID])*I;

                            const scalar sigmaEff = 
                            min(
                                sqrt((3.0/2.0)*magSqr(dev(sigmaPI))) 
                                - 3*muP[faceID]*DepsilonCreepEqP[faceID]
                                , 
                                3e10
                                );

                            const scalar DepsilonCreepEqPPrev = DepsilonCreepEqP[faceID];

                            //- Reference to temperature field [K]
                            const scalar T = Tp[faceID];

                            /*--- IRRADIATION CREEP RATE [1/hr] ---*/ 
                            //- Fast flux 
                            scalar phi = fastFluxP[faceID];

                            //- Tobbe correlation for 15-15 Ti austenitic [1/hr]
                            const scalar creepIrrRate = 
                            1e-2 * par1 * avgE * phi * sigmaEff/1e6 ;

                            DepsilonCreepIrrEqP[faceID] =
                            (1 - relax_)*DepsilonCreepIrrEqP[faceID] 
                            + relax_*creepIrrRate*deltaT;

                            /*--- THERMAL CREEP RATE [1/hr] ---*/ 

                            scalar creepThRate = 1e-2 *
                            par2 * exp(-par3/(par4*T)) * Foam::sinh(par5*sigmaEff/1e6/(par4*T));

                            DepsilonCreepThEqP[faceID] = 
                            (1 - relax_)*DepsilonCreepThEqP[faceID] 
                            + relax_*creepThRate*deltaT;

                            /*--- TOTAL CREEP RATE [1/hr] ---*/ 
                            const scalar fTrial = creepThRate + creepIrrRate;

                            DepsilonCreepEqP[faceID] = 
                            (1-relax_)*DepsilonCreepEqP[faceID]
                            + relax_*fTrial*deltaT;

                            epsilonCreepIrrEqP[faceID] = 
                            DepsilonCreepIrrEqP[faceID] 
                            + epsilonCreepIrrEqOldP[faceID] ;
                            
                            epsilonCreepIrrP[faceID] = 
                            epsilonCreepIrrOldP[faceID] 
                            + 1.5*DepsilonCreepIrrEqP[faceID]
                            *dev(sigmaPI - 2*muP[faceID]*DepsilonCreepP[faceID])/max(sigmaEff, SMALL);

                            epsilonCreepThEqP[faceID] = 
                            DepsilonCreepThEqP[faceID] 
                            + epsilonCreepThEqOldP[faceID] ;
                            
                            epsilonCreepThP[faceID] = 
                            epsilonCreepThOldP[faceID] 
                            + 1.5*DepsilonCreepThEqP[faceID]
                            *dev(sigmaPI - 2*muP[faceID]*DepsilonCreepP[faceID])/max(sigmaEff, SMALL);

                            epsilonCreepEqP[faceID] = 
                            DepsilonCreepEqP[faceID] + epsilonCreepEqOldP[faceID] ;

                            DepsilonCreepP[faceID] = 
                            (1 - relax_)*DepsilonCreepP[faceID]
                            + relax_*1.5*DepsilonCreepEqP[faceID]
                            *dev(sigmaPI - 2*muP[faceID]*DepsilonCreepP[faceID])/max(sigmaEff, SMALL) ;
                            
                            epsilonCreepP[faceID] = 
                            epsilonCreepOldP[faceID] + DepsilonCreepP[faceID];
                            
                            const scalar tolTr = 
                            min(
                                    1
                                    , 
                                    mag
                                    ( 
                                        (DepsilonCreepEqP[faceID] - DepsilonCreepEqPPrev)
                                        /(max(DepsilonCreepEqPPrev , SMALL))
                                    )
                                );

                            tol = min(tol, tolTr );
                        }
                    }
                }
            }

        }while(tol > 1);
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


