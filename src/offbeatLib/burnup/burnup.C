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

#include "burnup.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(burnup, 0);
    defineRunTimeSelectionTable(burnup, dictionary);
    addToRunTimeSelectionTable
    (
        burnup,
        burnup,
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::burnup::burnup
(
    const fvMesh& mesh,
    const materials& mat,
    const dictionary& burnupDict
)
:
    mesh_(mesh),
    mat_(mat),
    burnupDict_(burnupDict)
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::burnup>
Foam::burnup::New
(
    const fvMesh& mesh,
    const materials& mat,
    const dictionary& solverDict
)
{
    // Default pick type
    word type("none");

    // Prepare options dict
    dictionary burnupDict
    (
        solverDict.subOrEmptyDict("burnupOptions")
    );

    if (solverDict.found("burnup"))
    {
        if (solverDict.isDict("burnup"))
        {
            // Case 1: burnup is a dictionary
            burnupDict = solverDict.subDict("burnup");
            burnupDict.lookup("type") >> type;
        }
        else
        {
            // Case 2: burnup is a word (type)
            solverDict.lookup("burnup")  >> type;
            // burnupDict stays as burnupOptions
        }
    }

    // Map deprecated `fromLatestTime` into new name `constant`
    if (type == "fromLatestTime")
    {
        WarningIn("burnup::New(const fvMesh&, const materials&, const dictinoary&)")
            << "Type 'fromLatestTime' for `burnup` is deprecated. " 
            << "Please use 'constant' instead." << nl << endl;
        
        type = "constant";
    }
    
    Info << "Selecting burnup: " << type << endl;

    auto cstrIter
        = dictionaryConstructorTablePtr_->find(type);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("burnup::New(const fvMesh& mesh, volScalarField& Q)")
            << "Unknown burnup type "
            << type << nl << nl
            << "Valid types are:" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    if (debug)
    {
        Info<< "Selecting burnup type "
            << type << endl;
    }

    return autoPtr<burnup>(cstrIter()
        (
            mesh,
            mat,
            burnupDict
        ));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::burnup::~burnup()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// ************************************************************************* //
