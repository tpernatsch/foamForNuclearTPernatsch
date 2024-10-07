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

#include "heatCapacityModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(heatCapacityModel, 0);
    defineRunTimeSelectionTable(heatCapacityModel, dictionary);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatCapacityModel::heatCapacityModel
(
    const fvMesh& mesh,
    const dictionary& dict, 
    const word defaultModel
)
:
    mesh_(mesh)
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::heatCapacityModel>
Foam::heatCapacityModel::New
(
    const fvMesh& mesh,
    const dictionary& dict, 
    const word defaultModel  
)
{
    word heatCapacityModelName;

    heatCapacityModelName = 
    dict.lookupOrDefault<word>("heatCapacityModel", defaultModel);

    auto cstrIter
        = dictionaryConstructorTablePtr_->find(heatCapacityModelName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("heatCapacityModel::New(const fvMesh&, const dictionary&, const word)")
            << "Unknown heatCapacityModel type "
            << heatCapacityModelName << nl << nl
            << "Valid types are:" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    if (debug)
    {
        Info<< "Selecting heatCapacityModel type "
            << heatCapacityModelName << endl;
    }

    return autoPtr<heatCapacityModel>(cstrIter()(mesh, dict, defaultModel));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::heatCapacityModel::~heatCapacityModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// ************************************************************************* //
