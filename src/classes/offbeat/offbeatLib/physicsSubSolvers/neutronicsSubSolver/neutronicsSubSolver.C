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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::neutronicsSubSolver::neutronicsSubSolver
(
    const fvMesh& mesh,
    const burnup& bu,
    const dictionary& neutronicsOptDict
)
:     
    mesh_(mesh),
    burnup_(bu),
    neutronicsOptDict_(neutronicsOptDict),
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
    // Initialize type for neutronicsSubSolver class
    word type;

    dictionary neutronicsOptDict
    (
        solverDict.subOrEmptyDict("neutronicsSolverOption")
    );

    solverDict.lookup("neutronicsSolver") >> type;
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
            neutronicsOptDict
        ));
}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::neutronicsSubSolver::~neutronicsSubSolver()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// ************************************************************************* //

