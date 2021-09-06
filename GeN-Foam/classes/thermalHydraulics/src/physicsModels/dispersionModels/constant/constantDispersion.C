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

#include "FFPair.H"
#include "constantDispersion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dispersionModels
{
    defineTypeNameAndDebug(constant, 0);
    addToRunTimeSelectionTable
    (
        dispersionModel,
        constant,
        dispersionModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dispersionModels::constant::constant
(
    const FFPair& pair,
    const dictionary& dict,
    const objectRegistry& objReg
)
:
    dispersionModel
    (
        pair,
        dict,
        objReg
    ),
    dispersedPhaseName_(dict.get<word>("dispersedPhase")),
    value_
    (
        (pair_.fluid1().name() == dispersedPhaseName_) ? 
        dict.getOrDefault("value", 1.0) :
        1.0-dict.getOrDefault("value", 1.0)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::dispersionModels::constant::value(const label& celli) const
{
    return value_;
}

// ************************************************************************* //
