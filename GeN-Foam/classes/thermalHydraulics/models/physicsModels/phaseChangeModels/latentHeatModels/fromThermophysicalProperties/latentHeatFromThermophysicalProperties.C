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

#include "latentHeatFromThermophysicalProperties.H"
#include "fluid.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace latentHeatModels
{
    defineTypeNameAndDebug(latentHeatFromThermophysicalProperties, 0);
    addToRunTimeSelectionTable
    (
        latentHeatModel,
        latentHeatFromThermophysicalProperties,
        latentHeatModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::latentHeatModels::latentHeatFromThermophysicalProperties::
latentHeatFromThermophysicalProperties
(
    const dictionary& dict, 
    const fluid& fluid1,
    const fluid& fluid2,
    const volScalarField& p,
    const volScalarField& dmdt,
    const volScalarField& iT
)
:
    latentHeatModel
    (
        dict,
        fluid1,
        fluid2,
        p,
        dmdt, 
        iT
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::latentHeatModels::latentHeatFromThermophysicalProperties::correct()
{
    L_ = LSign_*(vapour_.thermo().hc() - liquid_.thermo().hc());

    this->adjust();
}

// ************************************************************************* //
