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
#include "ChenFlowEnhancementFactor.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace flowEnhancementFactorModels
{
    defineTypeNameAndDebug(Chen, 0);
    addToRunTimeSelectionTable
    (
        flowEnhancementFactorModel,
        Chen,
        flowEnhancementFactorModels
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::flowEnhancementFactorModels::Chen::Chen
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
    XLM_
    (
        pair.fluidRef().XLM()
    )
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::flowEnhancementFactorModels::Chen::value
(
    const label& celli
) const
{
    scalar invX(1.0/XLM_[celli]);
    return 
        (invX > 0.1002071798) ?
        min(2.35*pow(0.213+invX, 0.736), 50) :
        1.0;
}

// ************************************************************************* //
