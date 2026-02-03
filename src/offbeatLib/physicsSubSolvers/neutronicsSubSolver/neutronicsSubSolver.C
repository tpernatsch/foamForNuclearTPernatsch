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

#include "neutronicsSubSolver.H"
#include "fvm.H"
#include "fvc.H"
#include "calculatedFvPatchFields.H"
#include "PstreamReduceOps.H"
#include "Tuple2.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchField.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(neutronicsSubSolver, 0);
    defineRunTimeSelectionTable(neutronicsSubSolver, dictionary);   

    addToRunTimeSelectionTable
    (
        neutronicsSubSolver, 
        neutronicsSubSolver, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Protected Member Functions * * * * * * * * * * * * //

int Foam::neutronicsSubSolver::nCorrectors(const word fieldName)
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

Foam::neutronicsSubSolver::neutronicsSubSolver
(
    const fvMesh& mesh,
    const burnup& bu,
    const dictionary& neutronicsDict
)
:     
    mesh_(mesh),
    burnup_(bu),
    neutronicsDict_(neutronicsDict),
    residual_(1.0),
    initialResidual_(1.0),
    relResidual_(1.0)
{
}

Foam::autoPtr<Foam::neutronicsSubSolver>
Foam::neutronicsSubSolver::New
(
    const fvMesh& mesh,
    const burnup& bu,
    const dictionary& solverDict
)
{
    // Default pick type
    word type("none");

    // Prepare options dict
    dictionary neutronicsDict
    (
        solverDict.subOrEmptyDict("neutronicsSolverOption")
    );

    if (solverDict.found("neutronicsSolver"))
    {
        if (solverDict.isDict("neutronicsSolver"))
        {
            // Case 1: neutronicsSolver is a dictionary
            neutronicsDict = solverDict.subDict("neutronicsSolver");
            neutronicsDict.lookup("type") >> type;
        }
        else
        {
            // Case 2: neutronicsSolver is a word (type)
            solverDict.lookup("neutronicsSolver") >> type;
            // neutronicsDict stays as neutronicsSolverOptions
        }
    }
    
    // Map deprecated `fromLatestTime` into new name `none`
    if (type == "fromLatestTime")
    {
        WarningIn("thermalSubSolver::New(const fvMesh&, burnup&, const dictinoary&)")
            << "Type 'fromLatestTime' for `neutronicsSolver` is deprecated. " 
            << "Please use 'none' instead." << nl << endl;
        
        type = "none";
    }
    Info << "Selecting neutronicsSolver model " << type << endl;
    
    auto cstrIter
        = dictionaryConstructorTablePtr_->find(type);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("neutronicsSubSolver::New(const fvMesh& mesh, burnup& bu, const dictionary& solverDict)")
            << "Unknown neutronicsSubSolver type "
            << type << nl << nl
            << "Valid types are:" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    if (debug)
    {
        Info<< "Selecting neutronicsSubSolver type "
            << type << endl;
    }

    return autoPtr<neutronicsSubSolver>(cstrIter()
        (
            mesh, 
            bu, 
            neutronicsDict
        ));
}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::neutronicsSubSolver::~neutronicsSubSolver()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// ************************************************************************* //

