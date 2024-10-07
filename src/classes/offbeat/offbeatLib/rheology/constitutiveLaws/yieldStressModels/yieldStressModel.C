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

#include "yieldStressModel.H"
#include "addToRunTimeSelectionTable.H"
#include "globalFieldLists.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(yieldStressModel, 0);
    defineRunTimeSelectionTable(yieldStressModel, dictionary);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::yieldStressModel::yieldStressModel
(
    const fvMesh& mesh,
    const dictionary& lawDict
)
:
    mesh_(mesh),
    lawDict_(lawDict),
    sigmaY_(createOrLookup<scalar>(mesh, "sigmaY", dimPressure))
{
}

// * * * * * * * * * * * * * * * * Selectors  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::yieldStressModel>
Foam::yieldStressModel::New
(
    const fvMesh& mesh,
    const dictionary& lawDict
)
{
    word yieldStressModelName(lawDict.lookupOrDefault<word>("yieldStressModel", "hardening"));

    auto cstrIter
        = dictionaryConstructorTablePtr_->find(yieldStressModelName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("yieldStressModel::New(const fvMesh&, const dictionary&, const labelList&)")
            << "Unknown yieldStressModel type "
            << yieldStressModelName << nl << nl
            << "Valid types are:" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    if (debug)
    {
        Info<< "Selecting yieldStressModel type "
            << yieldStressModelName << endl;
    }

    return autoPtr<yieldStressModel>(cstrIter()(mesh, lawDict));
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::yieldStressModel::~yieldStressModel()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// ************************************************************************* //


