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

#include "FSPair.H"
#include "isomolarFluidDiameter.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fluidDiameterModels
{
    defineTypeNameAndDebug(isomolar, 0);
    addToRunTimeSelectionTable
    (
        fluidDiameterModel,
        isomolar,
        fluidDiameterModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidDiameterModels::isomolar::isomolar
(
    const FSPair& pair,
    const dictionary& dict,
    const objectRegistry& objReg
)
:
    fluidDiameterModel
    (
        pair,
        dict,
        objReg
    ),
    d0_(dict.get<scalar>("value")),
    p0_(dict.get<scalar>("p0")),
    T0_(dict.get<scalar>("T0")),
    p_(pair.mesh().lookupObject<volScalarField>("p")),
    T_(pair.fluidRef().thermo().T())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::fluidDiameterModels::isomolar::value
(
    const label& celli
) const
{
    return d0_*Foam::cbrt((p0_*T_[celli])/(p_[celli]*T0_));
}

// ************************************************************************* //
