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

#include "densificationModel.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchField.H"
#include "globalFieldLists.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(densificationModel, 0);
    defineRunTimeSelectionTable(densificationModel, dictionary);

    addToRunTimeSelectionTable
    (
        densificationModel, 
        densificationModel, 
        dictionary
    );

    const char* densificationModel::group_ = ("behavioralModels::densification");

    defineParameter(densificationModel, F_epsilonDensification_, 
        "F_epsilonDensification", (dimless), 1.0);

    defineParameter(densificationModel, delta_epsilonDensification_, 
        "delta_epsilonDensification", (dimless), 0.0);  
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::densificationModel::densificationModel
(
    const fvMesh& mesh,
    const dictionary& dict 
)
:
    mesh_(mesh),
    epsilonDensification_
    (
        createOrLookup<scalar>(mesh, "epsilonDensification", dimless, 0.0,
                zeroGradientFvPatchField<scalar>::typeName)
    )
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::densificationModel>
Foam::densificationModel::New
(
    const fvMesh& mesh,
    const dictionary& dict, 
    const word defaultModel        
)
{
    word densificationModelName;

    densificationModelName = 
    dict.lookupOrDefault<word>("densificationModel", defaultModel);

    auto cstrIter
        = dictionaryConstructorTablePtr_->find(densificationModelName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("densificationModel::New(const fvMesh&, const dictionary&, const word)")
            << "Unknown densificationModel type "
            << densificationModelName << nl << nl
            << "Valid types are:" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    if (debug)
    {
        Info<< "Selecting densificationModel type "
            << densificationModelName << endl;
    }

    return autoPtr<densificationModel>(cstrIter()(mesh, dict));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::densificationModel::~densificationModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// ************************************************************************* //
