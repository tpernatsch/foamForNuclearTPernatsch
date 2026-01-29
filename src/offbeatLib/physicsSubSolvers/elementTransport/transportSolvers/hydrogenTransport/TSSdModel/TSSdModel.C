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

#include "TSSdModel.H"
#include "addToRunTimeSelectionTable.H"
#include "globalFieldLists.H"
#include "materials.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(TSSdModel, 0);
    defineRunTimeSelectionTable(TSSdModel, dictionary);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::TSSdModel::TSSdModel
(
    const fvMesh& mesh,
    const materials& mat,
    const dictionary& lawDict
)
:
    mesh_(mesh),
    mat_(mat),
    lawDict_(lawDict),
    TSSd_(createOrLookup<scalar>(mesh, "TSSd", dimless))//,
    //TSSd0(0)
{
}

// * * * * * * * * * * * * * * * * Selectors  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::TSSdModel>
Foam::TSSdModel::New
(
    const fvMesh& mesh,
    const materials& mat,
    const dictionary& lawDict
)
{
    word TSSdModelName;
    TSSdModelName = lawDict.lookupOrDefault<word>("TSSdModel", "TSSdBison");

    auto cstrIter
        = dictionaryConstructorTablePtr_->find(TSSdModelName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("TSSdModel::New(const fvMesh&, const dictionary&, const labelList&)")
            << "Unknown TSSdModel type "
            << TSSdModelName << nl << nl
            << "Valid types are:" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    //Info << endl;
    Info << tab << tab << "Selecting TSSdModel type:  "
    << TSSdModelName << endl;

    return autoPtr<TSSdModel>(cstrIter()(mesh, mat, lawDict));
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::TSSdModel::~TSSdModel()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// ************************************************************************* //


