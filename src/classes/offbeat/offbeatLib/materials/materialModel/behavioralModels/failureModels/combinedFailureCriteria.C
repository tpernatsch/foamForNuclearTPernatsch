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

#include "combinedFailureCriteria.H"
#include "addToRunTimeSelectionTable.H"
#include "offbeatTime.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(combinedFailureCriteria, 0);
    addToRunTimeSelectionTable(failureModel, combinedFailureCriteria, dictionary);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::combinedFailureCriteria::combinedFailureCriteria
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    failureModel(mesh, dict)
{
    // Read list of words for combined models
    const List<word> criteria(failureModelDict_.lookup("criteria"));

    // Create failure model objects and put them in the list
    failureModelList_.clear();
    failureModelList_.resize(criteria.size());

    forAll(criteria, i)
    {
        failureModelList_.set
        (
            i,
            combinedFailureCriteria::New
            (
                mesh,
                dict,
                criteria[i]
            )
        );
    }
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::failureModel>
Foam::combinedFailureCriteria::New
(
    const fvMesh& mesh,
    const dictionary& dict, 
    const word defaultModel   
)
{
    // Select the name given as asrgument in the New method
    word failureModelName = defaultModel;

    if ( failureModelName == Foam::combinedFailureCriteria::typeName )
    {
        // FATAL ERROR : cannot select the typeName of this class
        FatalErrorIn("combinedFailureCriteria::New(const fvMesh&, const dictionary&, const word)")
        << "failureModel type " << failureModelName << " cannot be selected"
        << " as element of the criteria list for combined failure models."
        << exit(FatalError);
    }

    auto cstrIter
        = dictionaryConstructorTablePtr_->find(failureModelName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("failureModel::New(const fvMesh&, const dictionary&, const word)")
            << "Unknown failureModel type "
            << failureModelName << nl << nl
            << "Valid types are:" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    if (debug)
    {
        Info<< "Selecting failureModel type "
            << failureModelName << endl;
    }

    return autoPtr<failureModel>(cstrIter()(mesh, dict));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::combinedFailureCriteria::~combinedFailureCriteria()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::combinedFailureCriteria::isFailed
(
    const labelList& addr
)
{

    // Set failure switch to false
    failed_ = false;

    // Loop over all criteria
    forAll(failureModelList_, i)
    {
        failed_ = failureModelList_[i].isFailed(addr);
        
        if( failed_ )
        {
            return failed_;
        }
    }

    return failed_;

}


// ************************************************************************* //
