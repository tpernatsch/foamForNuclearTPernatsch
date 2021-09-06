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
#include "constantTemperatureSaturation.H"
#include "addToRunTimeSelectionTable.H"
#include "phaseChangeModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace saturationModels
{
    defineTypeNameAndDebug(constantTemperature, 0);
    addToRunTimeSelectionTable
    (
        saturationModel,
        constantTemperature,
        saturationModels
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::saturationModels::constantTemperature::
constantTemperature
(
    const phaseChangeModel& pcm,
    const dictionary& dict, 
    const objectRegistry& objReg
)
:
    saturationModel
    (
        pcm,
        dict,
        objReg
    ),
    TSat_(dict.get<scalar>("value")),
    p_(pcm.mesh().lookupObject<volScalarField>("p"))
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::saturationModels::constantTemperature::valueTSat
(
    const label& celli
) const
{
    return TSat_;
}

Foam::scalar Foam::saturationModels::constantTemperature::valuePSat
(
    const label& celli
) const
{
    return p_[celli];
}

Foam::scalar Foam::saturationModels::constantTemperature::valueLnPSat
(
    const label& celli
) const
{
    return log(p_[celli]);
}

Foam::scalar Foam::saturationModels::constantTemperature::valuePSatPrime
(
    const label& celli
) const
{
    return 0;
}

// ************************************************************************* //
