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

#include "fpYieldBurnupDependent.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(fpYieldBurnupDependent, 0);
    addToRunTimeSelectionTable(fpYieldModel, fpYieldBurnupDependent, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fpYieldBurnupDependent::fpYieldBurnupDependent
(
    const dictionary& dict,
    const scalar enrichment
)
:
    fpYieldModel(dict, enrichment),
    yieldCoeff_(readScalar(dict.lookup("yieldCoeff"))),
    burnupExponent_(readScalar(dict.lookup("burnupExponent"))),
    burnupConversion_(dict.lookupOrDefault<scalar>("burnupConversion", 8.232))
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fpYieldBurnupDependent::updateYield
(
    scalarField& Y,
    const scalarField& Bu
) const
{
    forAll(Y, cellI)
    {
        // Convert Bu [MWd/tUO2] to %FIMA, clamp to ≥ 1 % to avoid zero exponent
        const scalar b = max(1.0, Bu[cellI] * 1e-3 / burnupConversion_);
        Y[cellI] = yieldCoeff_ * pow(b, burnupExponent_);
    }
}

// ************************************************************************* //
