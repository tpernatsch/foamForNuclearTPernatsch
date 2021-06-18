/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "powerOffCriterionTimer.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace powerOffCriterionModels
{
    defineTypeNameAndDebug(timer, 0);
    addToRunTimeSelectionTable
    (
        powerOffCriterionModel, 
        timer, 
        powerOffCriterionModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::powerOffCriterionModels::timer::timer
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    powerOffCriterionModel
    (
        mesh,
        dict
    ),
    t0_(this->get<scalar>("time"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::powerOffCriterionModels::timer::~timer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::powerOffCriterionModels::timer::powerOffCriterion()
{   
    scalar t(mesh_.time().timeOutputValue());
    if (t >= t0_) 
    {
        Info << "Power off" << endl;
        return true;
    }
    return false;
}


// ************************************************************************* //
