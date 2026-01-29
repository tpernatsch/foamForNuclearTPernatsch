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

#include "PoissonRatioModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(PoissonRatioModel, 0);
    defineRunTimeSelectionTable(PoissonRatioModel, dictionary);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PoissonRatioModel::PoissonRatioModel
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

Foam::autoPtr<Foam::PoissonRatioModel>
Foam::PoissonRatioModel::New
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word defaultModel  
)
{
    // Default pick type
    word PoissonRatioModelName = defaultModel;

    // Prepare model dict
    dictionary modelDict
    (
        dict.subOrEmptyDict("PoissonRatio")
    );

    if(dict.found("PoissonRatioModel"))
    {
        dict.lookup("PoissonRatioModel") >> PoissonRatioModelName;

        // Revert modelDict to materialDict
        modelDict = dict;        
            
        WarningIn("PoissonRatioModel::New(const fvMesh&, const dictionary&, const word)")
            << "Keyword 'PoissonRatioModel' is deprecated. " 
            << "Please create a 'PoissonRatio' dictionary instead, where the typeName of the model is expressed via a keyword 'type'." << nl << endl;
    }
    else
    {
        if(dict.found("PoissonRatio"))
        {
            modelDict.lookup("type") >> PoissonRatioModelName;
        }
        else
        {
            // else the name remains 'defaultModel'
            modelDict = dict;
        }
    }

    auto cstrIter
        = dictionaryConstructorTablePtr_->find(PoissonRatioModelName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("PoissonRatioModel::New(const fvMesh&, const dictionary&, const word)")
            << "Unknown PoissonRatio type "
            << PoissonRatioModelName << nl << nl
            << "Valid types are:" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    if (debug)
    {
        Info<< "Selecting PoissonRatio type "
            << PoissonRatioModelName << endl;
    }

    return autoPtr<PoissonRatioModel>(cstrIter()(mesh, modelDict, dict, defaultModel));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PoissonRatioModel::~PoissonRatioModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// ************************************************************************* //
