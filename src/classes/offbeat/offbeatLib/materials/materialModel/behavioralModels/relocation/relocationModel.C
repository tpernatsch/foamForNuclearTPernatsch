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

#include "relocationModel.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchField.H"
#include "globalFieldLists.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(relocationModel, 0);
    defineRunTimeSelectionTable(relocationModel, dictionary);
    
    addToRunTimeSelectionTable
    (
        relocationModel, 
        relocationModel, 
        dictionary
    );

    const char* relocationModel::group_ = ("behavioralModels::relocation");

    defineParameter(relocationModel, F_epsilonRelocation_, 
        "F_epsilonRelocation", (dimless), 1.0);
    defineParameter(relocationModel, F_epsilonRelocationRecovery_, 
        "F_epsilonRelocationRecovery", (dimless), 1.0);

    defineParameter(relocationModel, delta_epsilonRelocation_, 
        "delta_epsilonRelocation", (dimless), 0.0);  
    defineParameter(relocationModel, delta_epsilonRelocationRecovery_, 
        "delta_epsilonRelocationRecovery", (dimless), 0.0);  
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relocationModel::relocationModel
(
    const fvMesh& mesh,
    const dictionary& dict  
)
:
    mesh_(mesh),
    epsilonRelocation_
    (
        createOrLookup<scalar>(mesh, "epsilonRelocation", dimless, 0.0,
                zeroGradientFvPatchField<scalar>::typeName)
    ),
    epsilonRecoveredRelocation_
    (
        createOrLookup<scalar>(mesh, "epsilonRecoveredRelocation", dimless, 0.0,
                zeroGradientFvPatchField<scalar>::typeName)
    )
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::relocationModel>
Foam::relocationModel::New
(
    const fvMesh& mesh,
    const dictionary& dict, 
    const word defaultModel        
)
{
    word relocationModelName;

    relocationModelName = 
    dict.lookupOrDefault<word>("relocationModel", defaultModel);

    auto cstrIter
        = dictionaryConstructorTablePtr_->find(relocationModelName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("relocationModel::New(const fvMesh&, const dictionary&, const word)")
            << "Unknown relocationModel type "
            << relocationModelName << nl << nl
            << "Valid types are:" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    if (debug)
    {
        Info<< "Selecting relocationModel type "
            << relocationModelName << endl;
    }

    return autoPtr<relocationModel>(cstrIter()(mesh, dict));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::relocationModel::~relocationModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// ************************************************************************* //
