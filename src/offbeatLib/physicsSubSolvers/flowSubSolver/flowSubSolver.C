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

#include "flowSubSolver.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(flowSubSolver, 0);
    defineRunTimeSelectionTable(flowSubSolver, dictionary);   

    addToRunTimeSelectionTable
    (
        flowSubSolver, 
        flowSubSolver, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::flowSubSolver::flowSubSolver
(
    const Time& runTime,
    const dictionary& dict
)
:     
    runTime_(runTime),
    flowSubSolverDict_(dict)
{}


Foam::autoPtr<Foam::flowSubSolver>
Foam::flowSubSolver::New
(
    const Time& runTime,
    const dictionary& solverDict
)
{
    // Read type
    word type = flowSubSolver::typeName;
    
    if (solverDict.found("flowSolver"))
    {
        solverDict.lookup("flowSolver") >> type;
        Info << "Selecting flow solver model " << type << endl;
    }
    else
    {
        Info << "No flow solver selected. Flow solver is disabled." << endl;
    }

    // Runtime lookup of type
    dictionaryConstructorTable::iterator cstrIter
        = dictionaryConstructorTablePtr_->find(type);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("flowSubSolver::New(fvMesh& mesh, ...)")
            << "Unknown flowSubSolver type "
            << type << nl << nl
            << "Valid types are:" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    // Construct and return object
    return autoPtr<flowSubSolver>(cstrIter()
        (
            runTime,
            solverDict.subOrEmptyDict("flowSolverOptions")
        ));
}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::flowSubSolver::~flowSubSolver()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// ************************************************************************* //

