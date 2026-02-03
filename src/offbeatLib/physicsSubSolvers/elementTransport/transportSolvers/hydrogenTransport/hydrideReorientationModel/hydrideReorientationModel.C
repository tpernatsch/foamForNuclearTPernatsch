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

#include "hydrideReorientationModel.H"
#include "addToRunTimeSelectionTable.H"
#include "globalFieldLists.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * Static Data Members  * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(hydrideReorientationModel, 0);
    defineRunTimeSelectionTable(hydrideReorientationModel, dictionary);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hydrideReorientationModel::hydrideReorientationModel
(
    const fvMesh& mesh,
    const dictionary& lawDict
)
:
    mesh_(mesh),
    lawDict_(lawDict),
    fracRHy_(createOrLookup<scalar>(mesh, "fracRHy", dimless, 0.0)),
    radialHydrides_(createOrLookup<scalar>(mesh, "radialHydrides", dimless, 0.0))
{}

// * * * * * * * * * * * * * * * * Selectors  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::hydrideReorientationModel>
Foam::hydrideReorientationModel::New
(
    const fvMesh& mesh,
    const dictionary& lawDict
)
{
    word hydrideReorientationModelName;

    hydrideReorientationModelName =
        lawDict.lookupOrDefault<word>("hydrideReorientationModel", "hydrideReorientationDesquines");

    auto cstrIter =
        dictionaryConstructorTablePtr_->find(hydrideReorientationModelName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("hydrideReorientationModel::New(const fvMesh&, const dictionary&, const labelList&)")
            << "Unknown hydrideReorientationModel type "
            << hydrideReorientationModelName << nl << nl
            << "Valid types are:" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    Info<< tab << "Selecting hydride reorientation model "
        << hydrideReorientationModelName << endl;

    return autoPtr<hydrideReorientationModel>(cstrIter()(mesh, lawDict));
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::hydrideReorientationModel::~hydrideReorientationModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// ************************************************************************ //
