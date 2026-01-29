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
    const dictionary& dict,
    const dictionary& baseDict 
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
    // Default pick type
    word densificationModelName = defaultModel;

    // Prepare model dict
    dictionary modelDict
    (
        dict.subOrEmptyDict("densification")
    );

    if(dict.found("densificationModel"))
    {
        dict.lookup("densificationModel") >> densificationModelName;

        // Revert modelDict to materialDict
        modelDict = dict;        
            
        WarningIn("densificationModel::New(const fvMesh&, const dictionary&, const word)")
            << "Keyword 'densificationModel' is deprecated. " 
            << "Please create a 'densification' dictionary instead, where the typeName of the model is expressed via a keyword 'type'." << nl << endl;
    }
    else
    {
        if(dict.found("densification"))
        {
            modelDict.lookup("type") >> densificationModelName;
        }
        else
        {
            // else the name remains 'defaultModel'
            modelDict = dict;
        }
    }

    auto cstrIter
        = dictionaryConstructorTablePtr_->find(densificationModelName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("densificationModel::New(const fvMesh&, const dictionary&, const word)")
            << "Unknown densification type "
            << densificationModelName << nl << nl
            << "Valid types are:" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    if (debug)
    {
        Info<< "Selecting densification type "
            << densificationModelName << endl;
    }

    return autoPtr<densificationModel>(cstrIter()(mesh, modelDict, dict));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::densificationModel::~densificationModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// ************************************************************************* //
