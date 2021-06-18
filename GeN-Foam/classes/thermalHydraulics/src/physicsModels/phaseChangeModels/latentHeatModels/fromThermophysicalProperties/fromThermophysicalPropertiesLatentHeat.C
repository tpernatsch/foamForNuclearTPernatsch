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

#include "FFPair.H"
#include "fromThermophysicalPropertiesLatentHeat.H"
#include "addToRunTimeSelectionTable.H"
#include "phaseChangeModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace latentHeatModels
{
    defineTypeNameAndDebug(fromThermophysicalPropertiesLatentHeat, 0);
    addToRunTimeSelectionTable
    (
        latentHeatModel,
        fromThermophysicalPropertiesLatentHeat,
        latentHeatModels
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::latentHeatModels::fromThermophysicalPropertiesLatentHeat::
fromThermophysicalPropertiesLatentHeat
(
    const phaseChangeModel& pcm,
    const dictionary& dict, 
    const objectRegistry& objReg
)
:
    latentHeatModel
    (
        pcm,
        dict,
        objReg
    )
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void 
Foam::latentHeatModels::fromThermophysicalPropertiesLatentHeat::correctField
(
    volScalarField& L
) const
{
    L = LSign_*(vapour_.thermo().hc() - liquid_.thermo().hc());
    this->adjust(L);
}

// ************************************************************************* //
