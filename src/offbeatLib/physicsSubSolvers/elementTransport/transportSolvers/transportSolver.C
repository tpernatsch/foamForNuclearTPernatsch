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
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

int Foam::transportSolver::nCorrectors(const word fieldName)
{
    const dictionary& stressControl =
        mesh_.solutionDict().subDict("stressAnalysis");

    int nCorrectors = 1;
    if(stressControl.found("nCorrectors"))
    {
        if(stressControl.isDict("nCorrectors"))
        {
        #ifdef OPENFOAMFOUNDATION
            if (stressControl.subDict("nCorrectors").found(fieldName))
            {
                nCorrectors = stressControl.subDict("nCorrectors").lookup<int>(fieldName);
            }
            else if (stressControl.subDict("nCorrectors").found("default"))
            {
                nCorrectors = stressControl.subDict("nCorrectors").lookup<int>("default");
            }
        #elif OPENFOAMESI
            if (stressControl.subDict("nCorrectors").found(fieldName))
            {
                nCorrectors = stressControl.subDict("nCorrectors").get<int>(fieldName);
            }
            else if (stressControl.subDict("nCorrectors").found("default"))
            {
                nCorrectors = stressControl.subDict("nCorrectors").get<int>("default");
            }
        #endif
            else
            {
                FatalErrorInFunction()
                << "Number of inner correctors missing for field \"" << fieldName << "\""
                << " in the fvSolution/stressAnalysis/nCorrectors sub-dict." << nl
                << "Please set this value using either the \"" << fieldName << "\"" 
                << " or \"default\" keywords" << nl << exit(FatalError);
            }
        }
        else
        {
        #ifdef OPENFOAMFOUNDATION
            nCorrectors = stressControl.lookup<int>("nCorrectors");
        #elif OPENFOAMESI
            nCorrectors = stressControl.get<int>("nCorrectors");
        #endif
        }
    }

    return nCorrectors;
}

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

    Info << tab << "Selecting transportSolver " << type << endl;
    
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

    bool useRelRes(stressControl.lookupOrDefault("useRelRes"+fieldName, true));
    
    scalar convergenceTolerance
    (
        stressControl.lookupOrDefault(fieldName, 1e-6)
    );

    if ( useRelRes )
    {
        scalar relConvergenceTolerance
        (
            stressControl.lookupOrDefault("rel"+fieldName, 1e-6)
        );

        // Consider converged if both initialResidual and relResidual are 
        // below threshold or if one of the two is 10 times below threshold.
        converged = 
        ( 
            initialResidual_ < convergenceTolerance 
            && 
            relResidual_ < relConvergenceTolerance 
        )
        or
        (
            initialResidual_ < convergenceTolerance/10
        )
        or
        (
            relResidual_ < relConvergenceTolerance/10 
        )
        ; 
    }
    else
    {
        // In this case only the residual computed by OpenFOAM is accounted for 
        // the convergence monitoring.
        converged = 
        ( 
            initialResidual_ < convergenceTolerance
        ); 
    }

    return converged;
}

// ************************************************************************* //

