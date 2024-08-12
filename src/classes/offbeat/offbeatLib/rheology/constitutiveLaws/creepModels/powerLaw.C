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

#include "powerLaw.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(powerLaw, 0);
    addToRunTimeSelectionTable
    (
        creepModel, 
        powerLaw, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::powerLaw::powerLaw
(
    const fvMesh& mesh,
    const dictionary& lawDict
)
:
    creepModel(mesh, lawDict),
    timeIndex_(0),
    relax_(lawDict.lookupOrDefault<scalar>("relax", 1.0 )),
    B(readScalar(lawDict.lookup("B"))),
    sigmaC(readScalar(lawDict.lookup("sigmaC"))),
    n(readScalar(lawDict.lookup("n")))
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::powerLaw::~powerLaw()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::powerLaw::correctCreep
(
    volSymmTensorField& epsilonEl,
    const labelList& addr
)
{
    const volScalarField& mu(mesh_.lookupObject<volScalarField>("mu"));
    const volScalarField& lambda(mesh_.lookupObject<volScalarField>("lambda"));
    const volScalarField& Temp(mesh_.lookupObject<volScalarField>("T"));
    
    //- Irradiation time step [hr]
    scalar deltaT(0);
    
    if(mesh_.time().timeIndex() > 0)
    {
        deltaT = mesh_.time().deltaTValue()/3600.0;

        const scalarField& epsilonCreepEqOld_ = 
        epsilonCreepEq_.oldTime().internalField();

        const symmTensorField& epsilonCreepOldI = 
        epsilonCreep_.oldTime().internalField();
    
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
                DepsilonCreepEq_[cellI] *= 0;

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

                            DepsilonCreepEqP[faceID] *= 0;
                        }
                    }
                }
            }

            timeIndex_ = mesh_.time().timeIndex();
        }

        scalar maxTemp(0);
        scalar maxStress(0);
        // scalar volTot(0);

        forAll(addr, addrI)
        {
            const label cellI = addr[addrI];

            const symmTensor sigmaI = 2*mu[cellI]*epsilonEl[cellI] + lambda[cellI]*tr(epsilonEl[cellI])*I;

            //- Calculate the von Mises stress [Pa]
            const scalar sigmaEff = 
            min
            (
                sqrt((3.0/2.0)
                *magSqr(dev(sigmaI))) 
                - 3*mu[cellI]*DepsilonCreepEq_[cellI] 
                , 3e10
            );
            
            maxStress = max(sigmaEff, maxStress);

            // const scalar DepsilonCreepEqPrev = DepsilonCreepEq_[cellI];

            //- Reference to temperature field [K]
            const scalar T = Temp.internalField()[cellI];

            maxTemp = max(T, maxTemp);

            // * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * #

            //- Calculate creep strain rate with power law
            const scalar creepRate = B * pow(sigmaEff*1e-6/sigmaC, n);

            //- Calculate creep strain
            DepsilonCreepEq_[cellI] = 
            (1-relax_)*DepsilonCreepEq_[cellI] + relax_*creepRate*deltaT;

            // * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * #

            epsilonCreepEq_[cellI] = 
            epsilonCreepEqOld_[cellI] + DepsilonCreepEq_[cellI];

            DepsilonCreep_[cellI] = 
            (1 - relax_)*DepsilonCreep_[cellI] + relax_*1.5*DepsilonCreepEq_[cellI]*dev(sigmaI 
                - 2*mu[cellI]*DepsilonCreep_[cellI])/max(sigmaEff, SMALL); 

            epsilonCreep_[cellI] = 
            epsilonCreepOldI[cellI] + DepsilonCreep_[cellI]; 

            const cell& c = mesh_.cells()[cellI];  

            forAll(c, faceI)
            {
                const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);

                {
                    if (patchID > -1 and epsilonEl.boundaryField()[patchID].size())
                    {
                        const label faceID = 
                        mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);

                        // const scalarField& Tp = Temp.boundaryField()[patchID];
                        const scalarField& muP = mu.boundaryField()[patchID];
                        const scalarField& lambdaP = lambda.boundaryField()[patchID];
                        const symmTensorField& epsilonElP = epsilonEl.boundaryField()[patchID];

                        symmTensorField& epsilonCreepP = epsilonCreep_.boundaryFieldRef()[patchID];
                        symmTensorField& DepsilonCreepP = DepsilonCreep_.boundaryFieldRef()[patchID];

                        scalarField& epsilonCreepEqP = epsilonCreepEq_.boundaryFieldRef()[patchID];

                        const scalarField& epsilonCreepEqOldP = epsilonCreepEq_.oldTime().boundaryField()[patchID];

                        const symmTensorField& epsilonCreepOldP = epsilonCreep_.oldTime().boundaryField()[patchID];

                        scalarField& DepsilonCreepEqP = DepsilonCreepEq_.boundaryFieldRef()[patchID];

                        //- Calculate the von Mises stress [Pa]
                        const symmTensor sigmaPI = 
                        2*muP[faceID]*epsilonElP[faceID] + lambdaP[faceID]*tr(epsilonElP[faceID])*I;

                        const scalar sigmaEff = min(sqrt((3.0/2.0)*magSqr(dev(sigmaPI))) - 3*muP[faceID]*DepsilonCreepEqP[faceID], 3e10);

                        // * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * #

                        //- Calculate creep strain rate with power law
                        const scalar creepRate = B * pow(sigmaEff*1e-6/sigmaC, n);

                        //- Calculate creep strain
                        DepsilonCreepEqP[faceID] = 
                        (1-relax_)*DepsilonCreepEqP[faceID] + relax_*creepRate*deltaT;

                        // * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * #

                        epsilonCreepEqP[faceID] = 
                        DepsilonCreepEqP[faceID] + epsilonCreepEqOldP[faceID] ;

                        DepsilonCreepP[faceID] = 
                        (1 - relax_)*DepsilonCreepP[faceID] + relax_*1.5*DepsilonCreepEqP[faceID]*dev(sigmaPI - 2*muP[faceID]*DepsilonCreepP[faceID])/max(sigmaEff, SMALL) ;
                        
                        epsilonCreepP[faceID] = 
                        epsilonCreepOldP[faceID] + DepsilonCreepP[faceID];
                    }
                }
            }
        }
    }    
    
    forAll(addr, addrI)
    {
        const label cellI = addr[addrI];

        epsilonEl[cellI] -= DepsilonCreep_[cellI];     
        
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


