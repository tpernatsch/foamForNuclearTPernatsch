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

#include "poreVelocityModel.H"
#include "addToRunTimeSelectionTable.H"
#include "globalFieldLists.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(poreVelocityModel, 0);
    defineRunTimeSelectionTable(poreVelocityModel, dictionary);
    addToRunTimeSelectionTable
    (
        poreVelocityModel, 
        poreVelocityModel, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::poreVelocityModel::poreVelocityModel
(
    const fvMesh& mesh, 
    const dictionary& dict,
    const word defaultModel
)
:
    mesh_(mesh),
    materialModelDict_(dict),
    poreVelocity_
    (
        createOrLookup<vector>
        (
            mesh, 
            "poreVelocity", 
            dimLength/dimTime, 
            vector::zero,
            "calculated"
        )
    ),
    alphaPoreVelocity_(dict.lookupOrDefault<scalar>("alphaPoreVelocity", 1.0))
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::poreVelocityModel>
Foam::poreVelocityModel::New
(
    const fvMesh& mesh, 
    const dictionary& dict,
    const word defaultModel
)
{
    
    const word poreVelocityModelName =
        dict.lookupOrDefault<word>("poreVelocityModel", defaultModel);

    auto cstrIter
        = dictionaryConstructorTablePtr_->find(poreVelocityModelName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("poreVelocityModel::New(const fvMesh&, const dictionary&)")
            << "Unknown poreVelocityModel type "
            << poreVelocityModelName << nl << nl
            << "Valid types are:" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    if (debug)
    {
        Info << "Selecting poreVelocityModel type "
             << poreVelocityModelName << endl;
    }

    autoPtr<Foam::poreVelocityModel> 
    poreVelocityModelPtr(cstrIter()(mesh, dict, defaultModel));

    return poreVelocityModelPtr;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::poreVelocityModel::~poreVelocityModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// ************************************************************************* //
