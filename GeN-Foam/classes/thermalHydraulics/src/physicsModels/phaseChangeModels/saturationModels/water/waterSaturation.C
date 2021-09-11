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
#include "waterSaturation.H"
#include "addToRunTimeSelectionTable.H"
#include "phaseChangeModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace saturationModels
{
    defineTypeNameAndDebug(water, 0);
    addToRunTimeSelectionTable
    (
        saturationModel,
        water,
        saturationModels
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::saturationModels::water::
water
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
    iT_(pcm.pair().iT()),
    p_(pcm.mesh().lookupObject<volScalarField>("p")),
    A_(1e6*2.590718143628726e-10),
    B_(273.159),
    C_(4.247368421052632)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::saturationModels::water::valuePSat
(
    const label& celli
) const
{
    const scalar& T(iT_[celli]);
    return A_*pow(T-B_, C_);
}

Foam::scalar Foam::saturationModels::water::valuePSatPrime
(
    const label& celli
) const
{
    const scalar& T(iT_[celli]);
    return (C_-1.0)*A_*pow(T-B_, C_-1.0);
}

Foam::scalar Foam::saturationModels::water::valueLnPSat
(
    const label& celli
) const
{
    return log(valuePSat(celli));
}

Foam::scalar Foam::saturationModels::water::valueTSat
(
    const label& celli
) const
{
    return pow(p_[celli]/A_, 1.0/C_) + B_;       
}

// ************************************************************************* //
