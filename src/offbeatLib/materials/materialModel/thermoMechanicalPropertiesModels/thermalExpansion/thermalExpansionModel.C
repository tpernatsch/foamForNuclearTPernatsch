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

#include "thermalExpansionModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(thermalExpansionModel, 0);
    defineRunTimeSelectionTable(thermalExpansionModel, dictionary);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalExpansionModel::thermalExpansionModel
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& baseDict,
    const word defaultModel
)
:
    mesh_(mesh),
    Tref_(293.0)
{
    if(dict.dictName() == "thermalExpansion")
    {
        if(dict.found("Tref"))
        {
            dict.lookup("Tref") >> Tref_;
        }
    }
    else
    {
        // For retrocompatibility
        if(dict.found("Tref"))
        {
            Tref_ = dimensionedScalar(dict.lookup("Tref")).value();
        }
    }
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::thermalExpansionModel>
Foam::thermalExpansionModel::New
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word defaultModel  
)
{
    // Default pick type
    word thermalExpansionModelName = defaultModel;

    // Prepare model dict
    dictionary modelDict
    (
        dict.subOrEmptyDict("thermalExpansion")
    );

    if(dict.found("thermalExpansionModel"))
    {
        dict.lookup("thermalExpansionModel") >> thermalExpansionModelName;

        // Revert modelDict to materialDict
        modelDict = dict;        
            
        WarningIn("thermalExpansionModel::New(const fvMesh&, const dictionary&, const word)")
            << "Keyword 'thermalExpansionModel' is deprecated. " 
            << "Please create a 'thermalExpansion' dictionary instead, where the typeName of the model is expressed via a keyword 'type'." << nl << endl;
    }
    else
    {
        if(dict.found("thermalExpansion"))
        {
            modelDict.lookup("type") >> thermalExpansionModelName;
        }
        else
        {
            // else the name remains 'defaultModel'
            modelDict = dict;
        }
    }

    auto cstrIter
        = dictionaryConstructorTablePtr_->find(thermalExpansionModelName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("thermalExpansionModel::New(const fvMesh&, const dictionary&, const word)")
            << "Unknown thermalExpansion type "
            << thermalExpansionModelName << nl << nl
            << "Valid types are:" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    if (debug)
    {
        Info<< "Selecting thermalExpansion type "
            << thermalExpansionModelName << endl;
    }

    return autoPtr<thermalExpansionModel>(cstrIter()(mesh, modelDict, dict, defaultModel));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::thermalExpansionModel::~thermalExpansionModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// ************************************************************************* //
