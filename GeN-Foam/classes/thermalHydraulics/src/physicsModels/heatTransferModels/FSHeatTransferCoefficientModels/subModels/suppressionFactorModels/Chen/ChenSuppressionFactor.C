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
#include "ChenSuppressionFactor.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace suppressionFactorModels
{
    defineTypeNameAndDebug(Chen, 0);
    addToRunTimeSelectionTable
    (
        suppressionFactorModel,
        Chen,
        suppressionFactorModels
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::suppressionFactorModels::Chen::Chen
(
    const FSPair& pair,
    const dictionary& dict,
    const objectRegistry& objReg
)
:
    suppressionFactorModel
    (
        pair,
        dict,
        objReg
    ),
    otherFluidPtr_(nullptr)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::suppressionFactorModels::Chen::value
(
    const label& celli
) const
{
    if (otherFluidPtr_ == nullptr)
    {
        HashTable<const fluid*> fluids(pair_.mesh().lookupClass<fluid>());
        otherFluidPtr_ = 
            (fluids[fluids.toc()[0]]->name() == pair_.fluidRef().name()) ?
            fluids[fluids.toc()[1]] : fluids[fluids.toc()[0]];
    }   

    const scalar& a(pair_.fluidRef().normalized()[celli]);
    const scalar& oa(otherFluidPtr_->normalized()[celli]);
    scalar Re2pi
    (
        max
        (
            mag
            (
                a*pair_.fluidRef().rho()[celli]*pair_.fluidRef().U()[celli]
            +   oa*otherFluidPtr_->rho()[celli]*otherFluidPtr_->U()[celli]
            )*pair_.structureRef().Dh()[celli]/
            (
                    a*pair_.fluidRef().mu()[celli]
                +   oa*otherFluidPtr_->mu()[celli]
            ),
            10.0
        )
    );

    return 1.0/(1.0+2.53e-06*pow(Re2pi, 1.17));
}

// ************************************************************************* //
