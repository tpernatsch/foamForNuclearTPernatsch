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

#include "fpYieldUO2.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(fpYieldUO2, 0);
    addToRunTimeSelectionTable(fpYieldModel, fpYieldUO2, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fpYieldUO2::fpYieldUO2
(
    const dictionary& dict,
    const scalar enrichment
)
:
    fpYieldModel(dict, enrichment),
    fpName_(dict.get<word>("fpName")),
    // Cs: enrichment-dependent (threshold 0.175)
    yieldCs_
    (
        enrichment < dict.lookupOrDefault<scalar>("enrichmentThreshold", 0.175)
        ? dict.lookupOrDefault<scalar>("yieldLEU", 0.135)
        : dict.lookupOrDefault<scalar>("yieldHEU", 0.16)
    ),
    // Ag: burnup-dependent power law
    yieldAgCoeff_(dict.lookupOrDefault<scalar>("yieldAgCoeff",    1.31625e-3)),
    yieldAgExp_  (dict.lookupOrDefault<scalar>("yieldAgExp",      0.55734)),
    yieldAgConv_ (dict.lookupOrDefault<scalar>("yieldAgConv",     8.232)),
    // Kr, Sr: constant
    yieldKr_(dict.lookupOrDefault<scalar>("yieldKr", 0.297)),
    yieldSr_(dict.lookupOrDefault<scalar>("yieldSr", 0.10))
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fpYieldUO2::updateYield
(
    scalarField& Y,
    const scalarField& Bu
) const
{
    if (fpName_ == "Cs")
    {
        Y = yieldCs_;
    }
    else if (fpName_ == "Ag")
    {
        forAll(Y, cellI)
        {
            const scalar b = max(1.0, Bu[cellI] * 1e-3 / yieldAgConv_);
            Y[cellI] = yieldAgCoeff_ * pow(b, yieldAgExp_);
        }
    }
    else if (fpName_ == "Kr")
    {
        Y = yieldKr_;
    }
    else if (fpName_ == "Sr")
    {
        Y = yieldSr_;
    }
    else
    {
        FatalErrorIn("Foam::fpYieldUO2::updateYield")
            << "No built-in UO2 yield for species " << fpName_
            << ". Use an explicit yield model (enrichmentDependent, burnupDependent, constant)."
            << exit(FatalError);
    }
}

// ************************************************************************* //
