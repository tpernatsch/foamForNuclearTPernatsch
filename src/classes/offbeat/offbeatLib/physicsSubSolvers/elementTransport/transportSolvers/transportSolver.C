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

#include "transportSolver.H"
#include "fvm.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(transportSolver, 0);
    defineRunTimeSelectionTable(transportSolver, dictionary);   

    addToRunTimeSelectionTable
    (
        transportSolver, 
        transportSolver, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::transportSolver::transportSolver
(
    fvMesh& mesh,
    const materials& mat,
    const dictionary& elementTransportDict,
    const word& solverName
)
:     
    mesh_(mesh),
    mat_(mat),
    dict_(elementTransportDict),
    residual_(1.0),
    initialResidual_(1.0),
    relResidual_(1.0)
{}

Foam::autoPtr<Foam::transportSolver>
Foam::transportSolver::New
(
    fvMesh& mesh,
    const materials& mat,
    const dictionary& elementTransportDict,
    const word& solverName
)
{
    // Initialize type 
    word type(solverName);

    Info << tab << "Selecting transportSolver model " << type;
    
    auto cstrIter
        = dictionaryConstructorTablePtr_->find(type);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("transportSolver::New(fvMesh& mesh, materials& mat, word& solverName)")
            << "Unknown transportSolver type "
            << type << nl << nl
            << "Valid types are:" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    if (debug)
    {
        Info << tab << "Selecting transportSolver type "
             << type << endl;
    }

    return autoPtr<transportSolver>(cstrIter()
        (
            mesh,
            mat,
            elementTransportDict,
            solverName
        ));
}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::transportSolver::~transportSolver()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::transportSolver::converged
(
    const word& fieldName
)
{
    bool converged(false);

    const dictionary& stressControl = 
        mesh_.solutionDict().subDict("stressAnalysis");

    scalar convergenceTolerance(stressControl.lookupOrDefault(fieldName, 1e-6));
    // scalar relConvergenceTolerance(stressControl.lookupOrDefault("rel"+fieldName, 1e-6));

    // Consider converge if initialResidual  is below threshold or 
    // if it is 10 times below threshold.
    converged = 
            ( 
                (initialResidual_ < convergenceTolerance) 
            ); 

    return converged;
}

// ************************************************************************* //

