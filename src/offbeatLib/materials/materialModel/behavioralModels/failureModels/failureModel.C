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

#include "failureModel.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchField.H"
#include "offbeatTime.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(failureModel, 0);
    defineRunTimeSelectionTable(failureModel, dictionary);

    addToRunTimeSelectionTable
    (
        failureModel, 
        failureModel, 
        dictionary
    );
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::failureModel::failureModel
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& baseDict  
)
:
    mesh_(mesh),
    matDict_(dict),
    failureModelDict_(dict.dictName() == "failure"?
        dict : dict.subOrEmptyDict("failureModelOptions")),
    failedMaterial_(nullptr),
    failed_(false),
    stopIfFailed_(failureModelDict_.lookupOrDefault<bool>("stopIfFailed", false)),
    mapper_(nullptr)
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::failureModel>
Foam::failureModel::New
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word defaultModel        
)
{
    // Default pick type
    word failureModelName = defaultModel;

    // Prepare model dict
    dictionary modelDict
    (
        dict.subOrEmptyDict("failure")
    );

    if(dict.found("failureModel"))
    {
        dict.lookup("failureModel") >> failureModelName;

        // Revert modelDict to materialDict
        modelDict = dict;        
            
        WarningIn("failureModel::New(const fvMesh&, const dictionary&, const word)")
            << "Keyword 'failureModel' is deprecated. " 
            << "Please create a 'failure' dictionary instead, where the typeName of the model is expressed via a keyword 'type'." << nl << endl;
    }
    else
    {
        if(dict.found("failure"))
        {
            modelDict.lookup("type") >> failureModelName;
        }
        else
        {
            // else the name remains 'defaultModel'
            modelDict = dict;
        }
    }

    auto cstrIter
        = dictionaryConstructorTablePtr_->find(failureModelName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("failureModel::New(const fvMesh&, const dictionary&, const word)")
            << "Unknown failure type "
            << failureModelName << nl << nl
            << "Valid types are:" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    if (debug)
    {
        Info<< "Selecting failure type "
            << failureModelName << endl;
    }

    return autoPtr<failureModel>(cstrIter()(mesh, modelDict, dict));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::failureModel::~failureModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::failureModel::initializeFailedMaterialField()
{
    failedMaterial_.set   
    (
        new volScalarField
        (
            IOobject
            (
                "failedMaterial",
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar("failedMaterial", dimless, 0),
            "calculated"
        )
    );
}

void Foam::failureModel::checkFailure
(
    const labelList& addr
)
{
    const offbeatTime& userRunTime(
        static_cast<const offbeatTime&>(mesh_.time()));

    // Call function to return bool to understand if material is failed or no
    // call it isFailed() returning a bool
    if ( isFailed(addr) )
    {
        if ( stopIfFailed_ )
        {
            // Force writing of results and end simulation
            const_cast<Time&>(mesh_.time()).writeAndEnd();
        }
        
        InfoIn("Foam::failureModel::checkFailure(const labelList& addr)") << nl  
            << tab << "Failure occurred in the material."            << nl
            << tab << "Failure criterion : " << criterionName()      << nl
            << tab << "Time of failure : "   << userRunTime.userTime() << userRunTime.unit() << nl;
    }
}

// ************************************************************************* //
