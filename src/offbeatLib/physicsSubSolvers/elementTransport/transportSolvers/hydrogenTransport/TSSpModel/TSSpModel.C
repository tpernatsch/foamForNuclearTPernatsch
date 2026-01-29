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

#include "TSSpModel.H"
#include "addToRunTimeSelectionTable.H"
#include "globalFieldLists.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(TSSpModel, 0);
    defineRunTimeSelectionTable(TSSpModel, dictionary);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::TSSpModel::TSSpModel
(
    const fvMesh& mesh,   
    const materials& mat,
    const dictionary& lawDict
)
:
    mesh_(mesh),
    mat_(mat),
    lawDict_(lawDict),
    TSSp_(createOrLookup<scalar>(mesh, "TSSp", dimless))
{
}

// * * * * * * * * * * * * * * * * Selectors  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::TSSpModel>
Foam::TSSpModel::New
(
    const fvMesh& mesh,
    const materials& mat,
    const dictionary& lawDict
)
{
    word TSSpModelName;
    TSSpModelName = lawDict.lookupOrDefault<word>("TSSpModel", "TSSpBison");

    auto cstrIter
        = dictionaryConstructorTablePtr_->find(TSSpModelName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("TSSpModel::New(const fvMesh&, const dictionary&, const labelList&)")
            << "Unknown TSSpModel type "
            << TSSpModelName << nl << nl
            << "Valid types are:" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    Info<< tab << tab << "Selecting TSSpModel type:  "
    << TSSpModelName << endl;

    return autoPtr<TSSpModel>(cstrIter()(mesh, mat, lawDict));
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::TSSpModel::~TSSpModel()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// ************************************************************************* //


