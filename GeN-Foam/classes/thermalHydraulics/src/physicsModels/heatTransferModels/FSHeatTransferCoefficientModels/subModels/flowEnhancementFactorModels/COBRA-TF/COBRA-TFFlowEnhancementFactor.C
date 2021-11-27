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
#include "COBRA-TFFlowEnhancementFactor.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace flowEnhancementFactorModels
{
    defineTypeNameAndDebug(COBRA_TF, 0);
    addToRunTimeSelectionTable
    (
        flowEnhancementFactorModel,
        COBRA_TF,
        flowEnhancementFactorModels
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::flowEnhancementFactorModels::COBRA_TF::COBRA_TF
(
    const FSPair& pair,
    const dictionary& dict,
    const objectRegistry& objReg
)
:
    flowEnhancementFactorModel
    (
        pair,
        dict,
        objReg
    ),
    otherFluidPtr_(nullptr)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::flowEnhancementFactorModels::COBRA_TF::value
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
    const scalar& Re1pi(pair_.Re()[celli]);
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

    return min(max(pow(Re2pi/Re1pi, 0.8), 1.0), maxValue_);
}

// ************************************************************************* //
