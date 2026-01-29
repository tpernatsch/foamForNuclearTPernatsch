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
    DepsilonCreepEq_(createOrLookup(mesh, "DepsilonCreepEq", dimless, 0.0))
{    
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
    // Default pick type
    word type;

    // Prepare model dict
    dictionary creepDict(lawDict);

    if(lawDict.found("creepModel"))
    {
        // Use deprecated names
        lawDict.lookup("creepModel") >> type;
        //creepDict stays lawDict       
            
        WarningIn("creep::New(const fvMesh&, const dictionary&)")
            << "Keywords 'creepModel' is deprecated. " 
            << "Please create a 'creep' dictionary instead, where the typeName of the model is expressed via a keyword 'type'." << nl << endl;
    }
    else
    {
        if(lawDict.found("creep"))
        {
            creepDict = lawDict.subDict("creep");
            creepDict.lookup("type") >> type;
        }
    }

    auto cstrIter
        = dictionaryConstructorTablePtr_->find(type);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("creepModel::New(const fvMesh&, const dictionary&, const labelList& )")
            << "Unknown creepModel type "
            << type << nl << nl
            << "Valid types are:" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    if (debug)
    {
        Info<< "Selecting creepModel type "
            << type << endl;
    }

    return autoPtr<creepModel>(cstrIter()(mesh, creepDict));
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
// ************************************************************************* //


