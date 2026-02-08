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

#include "fastFlux.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fastFlux, 0);
    defineRunTimeSelectionTable(fastFlux, dictionary);
    addToRunTimeSelectionTable
    (
        fastFlux, 
        fastFlux, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fastFlux::fastFlux
(
    const fvMesh& mesh,
    const materials& materials,
    const dictionary& fastFluxDict
)
:
    mesh_(mesh),
    mat_(materials),
    fastFluxDict_(fastFluxDict)
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fastFlux>
Foam::fastFlux::New
(
    const fvMesh& mesh,
    const materials& mat,
    const dictionary& solverDict
)
{
    // Default pick type
    word type("none");

    // Prepare options dict
    dictionary fastFluxDict
    (
        solverDict.subOrEmptyDict("fastFluxOptions")
    );

    if (solverDict.found("fastFlux"))
    {
        if (solverDict.isDict("fastFlux"))
        {
            // Case 1: fastFlux is a dictionary
            fastFluxDict = solverDict.subDict("fastFlux");
            fastFluxDict.lookup("type") >> type;
        }
        else
        {
            // Case 2: fastFlux is a word (type)
            solverDict.lookup("fastFlux") >> type;
            // fastFluxDict stays as fastFluxOptions
        }
    }

    // Map deprecated `fromLatestTime` into new name `constant`
    if (type == "fromLatestTime")
    {
        WarningIn("fastFlux::New(const fvMesh&, const materials&, const dictinoary&)")
            << "Type 'fromLatestTime' for `fastFlux` is deprecated. " 
            << "Please use 'constant' instead." << nl << endl;
        
        type = "constant";
    }

    Info << "Selecting fastFlux " << type << endl;

    auto cstrIter
        = dictionaryConstructorTablePtr_->find(type);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("fastFlux::New(const fvMesh& mesh)")
            << "Unknown fastFlux dependence type "
            << type << nl << nl
            << "Valid types are:" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    if (debug)
    {
        Info<< "Selecting fastFlux type "
            << type << endl;
    }

    return autoPtr<fastFlux>(cstrIter()
        (
            mesh, 
            mat, 
            fastFluxDict
        ));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fastFlux::~fastFlux()
{}


// ************************************************************************* //
