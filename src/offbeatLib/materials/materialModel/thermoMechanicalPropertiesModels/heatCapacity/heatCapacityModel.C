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

#include "heatCapacityModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(heatCapacityModel, 0);
    defineRunTimeSelectionTable(heatCapacityModel, dictionary);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatCapacityModel::heatCapacityModel
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

Foam::autoPtr<Foam::heatCapacityModel>
Foam::heatCapacityModel::New
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word defaultModel  
)
{
    // Default pick type
    word heatCapacityModelName = defaultModel;

    // Prepare model dict
    dictionary modelDict
    (
        dict.subOrEmptyDict("heatCapacity")
    );

    if(dict.found("heatCapacityModel"))
    {
        dict.lookup("heatCapacityModel") >> heatCapacityModelName;

        // Revert modelDict to materialDict
        modelDict = dict;        
            
        WarningIn("heatCapacityModel::New(const fvMesh&, const dictionary&, const word)")
            << "Keyword 'heatCapacityModel' is deprecated. " 
            << "Please create a 'heatCapacity' dictionary instead, where the typeName of the model is expressed via a keyword 'type'." << nl << endl;
    }
    else
    {
        if(dict.found("heatCapacity"))
        {
            modelDict.lookup("type") >> heatCapacityModelName;
        }
        else
        {
            // else the name remains 'defaultModel'
            modelDict = dict;
        }
    }

    auto cstrIter
        = dictionaryConstructorTablePtr_->find(heatCapacityModelName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("heatCapacityModel::New(const fvMesh&, const dictionary&, const word)")
            << "Unknown heatCapacity type "
            << heatCapacityModelName << nl << nl
            << "Valid types are:" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    if (debug)
    {
        Info<< "Selecting heatCapacity type "
            << heatCapacityModelName << endl;
    }

    return autoPtr<heatCapacityModel>(cstrIter()(mesh, modelDict, dict, defaultModel));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::heatCapacityModel::~heatCapacityModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// ************************************************************************* //
