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

#include "ChenFlowFactor.H"
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
    defineTypeNameAndDebug(Chen, 0);
    addToRunTimeSelectionTable
    (
        flowFactorModel, 
        Chen, 
        flowFactorModels
    );
}
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatTransferModels::nucleateBoilingSubModels::flowFactorModels::
Chen::Chen
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
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::heatTransferModels::nucleateBoilingSubModels::flowFactorModels::
Chen::correct()
{
    //- Ref to fluid Lockhart-Martinelli parameter
    const volScalarField& X(FSPair_.fluidRef().XLM());

    //- I changed the invX threshold from 0.1 to 0.1002071798 to avoid
    //  minimal discontinuities in F when invX crosses the threshold.
    //  I am limiting it to a max of 50 just in case
    forAll(cellList_, i)
    {
        const label& celli(cellList_[i]);
        scalar invX(1.0/X[celli]);
        flowFactor_[i] = 
            (invX > 0.1002071798) ?
            min(2.35*pow(0.213+invX, 0.736), 50) :
            1.0;
    }
}


// ************************************************************************* //
