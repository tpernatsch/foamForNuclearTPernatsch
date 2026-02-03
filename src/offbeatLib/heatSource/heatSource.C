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

#include "heatSource.H"
#include "zeroGradientFvPatchField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(heatSource, 0);
    defineRunTimeSelectionTable(heatSource, dictionary);
    addToRunTimeSelectionTable
    (
        heatSource, 
        heatSource, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatSource::heatSource
(
    const fvMesh& mesh,
    const materials& materials,
    const dictionary& heatSourceDict
)
:
    mesh_(mesh),
    mat_(materials),
    heatSourceDict_(heatSourceDict)
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::heatSource>
Foam::heatSource::New
(
    const fvMesh& mesh,
    const materials& mat,
    const dictionary& solverDict
)
{
    // Default pick type
    word type("none");

    // Prepare options dict
    dictionary heatSourceDict
    (
        solverDict.subOrEmptyDict("heatSourceOptions")
    );

    if (solverDict.found("heatSource"))
    {
        if (solverDict.isDict("heatSource"))
        {
            // Case 1: heatSource is a dictionary
            heatSourceDict = solverDict.subDict("heatSource");
            heatSourceDict.lookup("type") >> type;
        }
        else
        {
            // Case 2: heatSource is a word (type)
            solverDict.lookup("heatSource") >> type;
            // heatSourceDict stays as heatSourceOptions
        }
    }

    // Map deprecated `fromLatestTime` into new name `constant`
    if (type == "fromLatestTime")
    {
        WarningIn("burnup::New(const fvMesh&, const materials&, const dictinoary&)")
            << "Type 'fromLatestTime' for `heatSource` is deprecated. " 
            << "Please use 'constant' instead." << nl << endl;
        
        type = "constant";
    }
    Info << "Selecting heatSource " << type << endl;

    auto cstrIter
        = dictionaryConstructorTablePtr_->find(type);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("heatSource::New(const fvMesh& mesh)")
            << "Unknown heatSource dependence type "
            << type << nl << nl
            << "Valid types are:" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    if (debug)
    {
        Info<< "Selecting heatSource type "
            << type << endl;
    }
    
    return autoPtr<heatSource>(cstrIter()
        (
            mesh, 
            mat, 
            heatSourceDict
        ));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::heatSource::~heatSource()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// ************************************************************************* //
