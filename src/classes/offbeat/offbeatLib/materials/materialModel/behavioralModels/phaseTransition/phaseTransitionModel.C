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

#include "phaseTransitionModel.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchField.H"
#include "calculatedFvPatchField.H"
#include "globalFieldLists.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseTransitionModel, 0);
    defineRunTimeSelectionTable(phaseTransitionModel, dictionary);
    
    addToRunTimeSelectionTable
    (
        phaseTransitionModel, 
        phaseTransitionModel, 
        dictionary
    );

    const char* phaseTransitionModel::group_ = ("behavioralModels::phaseTransition");

    defineParameter(phaseTransitionModel, F_betaFraction_, 
        "F_betaFraction", (dimless), 1.0);

    defineParameter(phaseTransitionModel, delta_betaFraction_, 
        "delta_betaFraction", (dimless), 0.0);  
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseTransitionModel::phaseTransitionModel
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    betaFraction_
    (
        createOrLookup<scalar>(mesh, "betaFraction", dimless, 0.0,
                calculatedFvPatchField<scalar>::typeName)
    )
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::phaseTransitionModel>
Foam::phaseTransitionModel::New
(
    const fvMesh& mesh,
    const dictionary& dict, 
    const word defaultModel        
)
{
    word phaseTransitionModelName;

    phaseTransitionModelName = 
    dict.lookupOrDefault<word>("phaseTransitionModel", defaultModel);

    auto cstrIter
        = dictionaryConstructorTablePtr_->find(phaseTransitionModelName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("phaseTransitionModel::New(const fvMesh&, const dictionary&, const word)")
            << "Unknown phaseTransitionModel type "
            << phaseTransitionModelName << nl << nl
            << "Valid types are:" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    if (debug)
    {
        Info<< "Selecting phaseTransitionModel type "
            << phaseTransitionModelName << endl;
    }

    return autoPtr<phaseTransitionModel>(cstrIter()(mesh, dict));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseTransitionModel::~phaseTransitionModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// ************************************************************************* //
