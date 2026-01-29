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

#include "densityModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(densityModel, 0);
    defineRunTimeSelectionTable(densityModel, dictionary);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::densityModel::densityModel
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& baseDict,
    const word defaultModel
)
:
    mesh_(mesh)
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::densityModel>
Foam::densityModel::New
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word defaultModel  
)
{
    // Default pick type
    word densityModelName = defaultModel;

    // Prepare model dict
    dictionary modelDict
    (
        dict.subOrEmptyDict("density")
    );

    if(dict.found("densityModel"))
    {
        dict.lookup("densityModel") >> densityModelName;

        // Revert modelDict to materialDict
        modelDict = dict;        
            
        WarningIn("densityModel::New(const fvMesh&, const dictionary&, const word)")
            << "Keyword 'densityModel' is deprecated. " 
            << "Please create a 'density' dictionary instead, where the typeName of the model is expressed via a keyword 'type'." << nl << endl;
    }
    else
    {
        if(dict.found("density"))
        {
            modelDict.lookup("type") >> densityModelName;
        }
        else
        {
            // else the name remains 'defaultModel'
            modelDict = dict;
        }
    }

    auto cstrIter
        = dictionaryConstructorTablePtr_->find(densityModelName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("densityModel::New(const fvMesh&, const dictionary&, const word)")
            << "Unknown density type "
            << densityModelName << nl << nl
            << "Valid types are:" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    if (debug)
    {
        Info<< "Selecting density type "
            << densityModelName << endl;
    }

    return autoPtr<densityModel>(cstrIter()(mesh, modelDict, dict, defaultModel));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::densityModel::~densityModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// ************************************************************************* //
