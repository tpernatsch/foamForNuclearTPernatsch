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

#include "thermalExpansionModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(thermalExpansionModel, 0);
    defineRunTimeSelectionTable(thermalExpansionModel, dictionary);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalExpansionModel::thermalExpansionModel
(
    const fvMesh& mesh,
    const dictionary& dict, 
    const word defaultModel
)
:
    mesh_(mesh),
    Tref_(dict.lookup("Tref"))
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::thermalExpansionModel>
Foam::thermalExpansionModel::New
(
    const fvMesh& mesh,
    const dictionary& dict, 
    const word defaultModel  
)
{
    word thermalExpansionModelName;

    thermalExpansionModelName = 
    dict.lookupOrDefault<word>("thermalExpansionModel", defaultModel);


    auto cstrIter
        = dictionaryConstructorTablePtr_->find(thermalExpansionModelName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("thermalExpansionModel::New(const fvMesh&, const dictionary&, const word)")
            << "Unknown thermalExpansionModel type "
            << thermalExpansionModelName << nl << nl
            << "Valid types are:" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    if (debug)
    {
        Info<< "Selecting thermalExpansionModel type "
            << thermalExpansionModelName << endl;
    }

    return autoPtr<thermalExpansionModel>(cstrIter()(mesh, dict, defaultModel));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::thermalExpansionModel::~thermalExpansionModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// ************************************************************************* //
