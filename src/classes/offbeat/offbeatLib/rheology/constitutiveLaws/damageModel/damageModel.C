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

#include "damageModel.H"
#include "addToRunTimeSelectionTable.H"
#include "globalFieldLists.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(damageModel, 0);
    defineRunTimeSelectionTable(damageModel, dictionary);

    addToRunTimeSelectionTable
    (
        damageModel, 
        damageModel,
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::damageModel::damageModel
(
    const fvMesh& mesh,
    const dictionary& lawDict
)
:
    mesh_(mesh),
    lawDict_(lawDict)
{
}

// * * * * * * * * * * * * * * * * Selectors  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::damageModel>
Foam::damageModel::New
(
    const fvMesh& mesh,
    const dictionary& lawDict
)
{
    word damageModelName(lawDict.lookupOrDefault<word>("damageModel", "none"));

   auto cstrIter
        = dictionaryConstructorTablePtr_->find(damageModelName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("damageModel::New(const fvMesh&, const dictionary&, const labelList&)")
            << "Unknown damageModel type "
            << damageModelName << nl << nl
            << "Valid types are:" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    if (debug)
    {
        Info<< "Selecting damageModel type "
            << damageModelName << endl;
    }

    return autoPtr<damageModel>(cstrIter()(mesh, lawDict));
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::damageModel::~damageModel()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// ************************************************************************* //


