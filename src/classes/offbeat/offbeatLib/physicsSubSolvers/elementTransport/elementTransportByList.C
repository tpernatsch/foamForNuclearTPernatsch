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

#include "elementTransportByList.H"
#include "zeroGradientFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "emptyFvPatchFields.H"
#include "gapGasModel.H"

#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(elementTransportByList, 0);
    addToRunTimeSelectionTable
    (
        elementTransport, 
        elementTransportByList, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::elementTransportByList::elementTransportByList
(
    fvMesh& mesh,
    const materials& mat,
    const dictionary& elementTransportDict
)
:
    elementTransport(mesh, mat, elementTransportDict),
    solverList_()
{

    const List<word> solvers(elementTransportDict.lookup("solvers"));

    solverList_.clear();
    solverList_.resize(solvers.size());

    forAll(solvers, i)
    {
        const word solver(solvers[i]);

        solverList_.set
        (
            i,
            transportSolver::New
            (
                mesh,
                mat,
                elementTransportDict,
                solver
            )
        );
        
        Info << endl;
    }
    
    Info << endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::elementTransportByList::~elementTransportByList()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::elementTransportByList::correct()
{

    forAll(solverList_, i)
    {
        solverList_[i].correct();
    }
}

bool Foam::elementTransportByList::converged()
{
    forAll(solverList_, i)
    {
        if( !solverList_[i].converged() )
        {
            return false;
        }   
    }

    return true;
}

Foam::scalar Foam::elementTransportByList::nextDeltaT()
{
    scalar nextDeltaT(GREAT);
    
    forAll(solverList_, i)
    {
        nextDeltaT = min
        (
            solverList_[i].nextDeltaT(),
            nextDeltaT
        );
    }


    return nextDeltaT;
}

void Foam::elementTransportByList::updateMesh()
{

    forAll(solverList_, i)
    {
        solverList_[i].updateMesh();
    }
}

// ************************************************************************* //
