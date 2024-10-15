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

#include "swellingModel.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchField.H"
#include "globalFieldLists.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(swellingModel, 0);
    defineRunTimeSelectionTable(swellingModel, dictionary);
    
    addToRunTimeSelectionTable
    (
        swellingModel, 
        swellingModel, 
        dictionary
    );

    const char* swellingModel::group_ = ("behavioralModels::swelling");

    defineParameter(swellingModel, F_epsilonSwelling_, 
        "F_epsilonSwelling", (dimless), 1.0);

    defineParameter(swellingModel, delta_epsilonSwelling_, 
        "delta_epsilonSwelling", (dimless), 0.0);  
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::swellingModel::swellingModel
(
    const fvMesh& mesh,
    const dictionary& dict   
)
:
    mesh_(mesh),
    epsilonSwelling_
    (
        createOrLookup<symmTensor>(mesh, "epsilonSwelling", dimless, symmTensor::zero,
                zeroGradientFvPatchField<symmTensor>::typeName)
    )
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::swellingModel>
Foam::swellingModel::New
(
    const fvMesh& mesh,
    const dictionary& dict, 
    const word defaultModel        
)
{
    word swellingModelName;

    swellingModelName = 
    dict.lookupOrDefault<word>("swellingModel", defaultModel);

    auto cstrIter
        = dictionaryConstructorTablePtr_->find(swellingModelName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("swellingModel::New(const fvMesh&, const dictionary&, const word)")
            << "Unknown swellingModel type "
            << swellingModelName << nl << nl
            << "Valid types are:" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    if (debug)
    {
        Info<< "Selecting swellingModel type "
            << swellingModelName << endl;
    }

    return autoPtr<swellingModel>(cstrIter()(mesh, dict));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::swellingModel::~swellingModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// ************************************************************************* //
