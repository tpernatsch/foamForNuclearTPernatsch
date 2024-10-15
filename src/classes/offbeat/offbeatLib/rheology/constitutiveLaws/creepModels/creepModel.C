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

#include "creepModel.H"
#include "addToRunTimeSelectionTable.H"
#include "globalFieldLists.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(creepModel, 0);
    defineRunTimeSelectionTable(creepModel, dictionary);

    addToRunTimeSelectionTable
    (
        creepModel, 
        creepModel,
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::creepModel::creepModel
(
    const fvMesh& mesh,
    const dictionary& lawDict
)
:
    mesh_(mesh),
    lawDict_(lawDict),
    epsilonCreep_(createOrLookup(mesh, "epsilonCreep", dimless, symmTensor::zero)),
    DepsilonCreep_(createOrLookup(mesh, "DepsilonCreep", dimless, symmTensor::zero)),
    epsilonCreepEq_(createOrLookup(mesh, "epsilonCreepEq", dimless, 0.0)),
    DepsilonCreepEq_(createOrLookup(mesh, "DepsilonCreepEq", dimless, 0.0)),
    maxAverageCreep_(GREAT),
    maxMaximumCreep_(GREAT)
{
    const dictionary& controlDict(mesh_.time().controlDict());

    bool adjustTime = controlDict.lookupOrDefault<bool>("adjustableTimeStep", false);

    if(adjustTime)
    {
        maxAverageCreep_ = controlDict.lookupOrDefault<scalar>(
            "maxAverageCreep", GREAT);
        maxMaximumCreep_ = controlDict.lookupOrDefault<scalar>(
            "maxMaximumCreep", GREAT);
    }
    
    epsilonCreep_.oldTime();
    DepsilonCreep_.oldTime();
    epsilonCreepEq_.oldTime();
    DepsilonCreepEq_.oldTime();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::creepModel::~creepModel()
{}


// * * * * * * * * * * * * * * * * Selectors  * * * * * * * * * * * * * * * * //


Foam::autoPtr<Foam::creepModel>
Foam::creepModel::New
(
    const fvMesh& mesh,
    const dictionary& lawDict
)
{
    word lawName;
    lawDict.lookup("creepModel") >> lawName;

    auto cstrIter
        = dictionaryConstructorTablePtr_->find(lawName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("creepModel::New(const fvMesh&, const dictionary&, const labelList& )")
            << "Unknown creepModel type "
            << lawName << nl << nl
            << "Valid types are:" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    if (debug)
    {
        Info<< "Selecting creepModel type "
            << lawName << endl;
    }

    return autoPtr<creepModel>(cstrIter()(mesh, lawDict));
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::creepModel::correctCreep
(
    volSymmTensorField& epsilonEl,
    const labelList& addr
)
{
    if(mesh_.time().timeIndex() > 0)
    {
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
                const label patchID(mesh_.boundaryMesh().whichPatch(c[faceI]));

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
}

void Foam::creepModel::correctAxialStrain
(
    volSymmTensorField& epsilon,
    const labelList& addr
)
{
    forAll(addr, i)
    {
        const label cellI = addr[i];

        epsilon[cellI].replace
        (
            symmTensor::ZZ,
            epsilon[cellI].component(symmTensor::ZZ)
           -(
                epsilonCreep_.internalField()[cellI].component(symmTensor::XX) + 
                epsilonCreep_.internalField()[cellI].component(symmTensor::YY)
            )
        );

        const cell& c = mesh_.cells()[cellI];  

        forAll(c, faceI)
        {
            const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);

            {
                if (patchID > -1 and epsilon.boundaryField()[patchID].size())
                {
                    const label faceID = 
                    mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);

                    symmTensorField& epsilonP = 
                    epsilon.boundaryFieldRef()[patchID];

                    const symmTensorField& epsilonCreepP = 
                    epsilonCreep_.boundaryField()[patchID];

                    epsilonP[faceID].replace
                    (
                        symmTensor::ZZ,
                        epsilonP[faceID].component(symmTensor::ZZ)
                       -(
                            epsilonCreepP[faceID].component(symmTensor::XX) + 
                            epsilonCreepP[faceID].component(symmTensor::YY)
                        )
                    );
                }
            }
        }
    }
}


Foam::scalar Foam::creepModel::nextDeltaT
(
    const labelList& addr
)
{
    //- Store old value of deltaT
    scalar oldDeltaT = mesh_.time().deltaT().value();

    //- Prepare new deltaT
    // scalar deltaT(0);

    //- Prepare average and maximum DepsilonCreepEq from last time step
    scalar averageDepsilonCreepEq(0);
    scalar maxDepsilonCreepEq(0);
    scalar volTot(0);

    const scalarField& DepsilonCreepEqI(DepsilonCreepEq_());
    const scalarField& V(mesh_.V());

    //- Find maximum and average DepsilonCreepEq
    forAll(addr, addrI)
    {
        const label cellI = addr[addrI];
        averageDepsilonCreepEq += DepsilonCreepEqI[cellI]*V[cellI];
        volTot += V[cellI];

        maxDepsilonCreepEq = max(maxDepsilonCreepEq, DepsilonCreepEqI[cellI]);
    }

    reduce(volTot, sumOp<scalar>());

    reduce(averageDepsilonCreepEq, sumOp<scalar>());

    reduce(maxDepsilonCreepEq, maxOp<scalar>());

    averageDepsilonCreepEq /= volTot;

    //- Calculate average and maximum creep rate from last time step
    scalar averageCreepRate(averageDepsilonCreepEq/oldDeltaT);
    scalar maxCreepRate(maxDepsilonCreepEq/oldDeltaT);

    Info << "deltaT: averageDepsilonCreepEq " << averageDepsilonCreepEq << endl;
    Info << "deltaT: averageCreepRate " << averageCreepRate << endl;
    Info << "deltaT: maxDepsilonCreepEq " << maxDepsilonCreepEq << endl;
    Info << "deltaT: maxCreepRate " << maxCreepRate << endl;

    return min
    (
        maxAverageCreep_/max(averageCreepRate, SMALL), 
        maxMaximumCreep_/max(maxCreepRate, SMALL)
    );


}
// ************************************************************************* //


