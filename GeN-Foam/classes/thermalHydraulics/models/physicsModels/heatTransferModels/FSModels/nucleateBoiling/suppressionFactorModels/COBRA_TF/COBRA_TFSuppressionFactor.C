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

#include "COBRA_TFSuppressionFactor.H"
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
namespace suppressionFactorModels
{
    defineTypeNameAndDebug(COBRA_TF, 0);
    addToRunTimeSelectionTable
    (
        suppressionFactorModel, 
        COBRA_TF, 
        suppressionFactorModels
    );
}
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatTransferModels::nucleateBoilingSubModels::suppressionFactorModels::
COBRA_TF::COBRA_TF
(
    const heatTransferModel& htm,
    const objectRegistry& objReg,
    const FSPair& FSPair
)
:
    suppressionFactorModel
    (
        htm,
        objReg,
        FSPair
    ),
    Twall_(FSPair.structure().Twall()),
    Tf_(FSPair.fluidRef().thermo().T()),
    Tsat_
    (
        //- The interface is always at saturation when doing simulations
        //  with phase change enabled, so use this for Tsat_
        FSPair.fluidRef().mesh().lookupObject<volScalarField>("T.interface")
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::heatTransferModels::nucleateBoilingSubModels::suppressionFactorModels::
COBRA_TF::correct()
{
    forAll(cellList_, i)
    {
        const label& celli(cellList_[i]);
        const scalar& Tf(Tf_[celli]);
        suppressionFactor_[i] = 
            max(Tf-Tsat_[celli], 0)/max(Twall_[celli]-Tf, 1e-3);
    }
}


// ************************************************************************* //
