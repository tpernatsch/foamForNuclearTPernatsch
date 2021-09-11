/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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

#include "forcedConstantPhaseChange.H"
#include "FFPair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseChangeModels
{
    defineTypeNameAndDebug(forcedConstantPhaseChange, 0);
    addToRunTimeSelectionTable
    (
        phaseChangeModel,
        forcedConstantPhaseChange,
        phaseChangeModels
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeModels::forcedConstantPhaseChange::forcedConstantPhaseChange
(
    FFPair& pair,
    const dictionary& dict
)
:
    phaseChangeModel
    (
        pair,
        dict
    )
{
    scalar value(this->get<scalar>("value"));
    wordList regions(this->get<wordList>("regions"));
    forAll(regions, i)
    {
        const labelList& zoneCellList(pair.mesh().cellZones()[regions[i]]);
        forAll(zoneCellList, j)
        {
            const label& celli(zoneCellList[j]);
            dmdt_[celli] = value;
        }
    }
    dmdt_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void 
Foam::phaseChangeModels::forcedConstantPhaseChange::correctInterfacialDmdt() 
{
}

// ************************************************************************* //
