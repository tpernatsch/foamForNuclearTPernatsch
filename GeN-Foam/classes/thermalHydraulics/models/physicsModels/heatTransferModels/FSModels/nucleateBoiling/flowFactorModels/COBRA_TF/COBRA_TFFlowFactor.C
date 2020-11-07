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

#include "COBRA_TFFlowFactor.H"
#include "addToRunTimeSelectionTable.H"
#include "fluid.H"
#include "structureModel.H"
#include "FSPair.H"
#include "heatTransferModel.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferModels
{
namespace nucleateBoilingSubModels
{
namespace flowFactorModels
{
    defineTypeNameAndDebug(COBRA_TF, 0);
    addToRunTimeSelectionTable
    (
        flowFactorModel, 
        COBRA_TF, 
        flowFactorModels
    );
}
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatTransferModels::nucleateBoilingSubModels::flowFactorModels::
COBRA_TF::COBRA_TF
(
    const heatTransferModel& htm,
    const objectRegistry& objReg,
    const FSPair& FSPair
)
:
    flowFactorModel
    (
        htm,
        objReg,
        FSPair
    ),
    Re1p_(FSPair.Re()),
    Re2p_(FSPair.fluidRef().mesh().lookupObject<volScalarField>("ReTwoPhase"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::heatTransferModels::nucleateBoilingSubModels::flowFactorModels::
COBRA_TF::correct()
{
    forAll(cellList_, i)
    {
        const label& celli(cellList_[i]);
        flowFactor_[i] = 
            max
            (
                pow(Re2p_[celli]/Re1p_[celli], 0.8),
                1.0
            );
    }
}


// ************************************************************************* //
