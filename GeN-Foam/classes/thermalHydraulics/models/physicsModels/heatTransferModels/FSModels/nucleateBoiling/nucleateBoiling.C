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

#include "nucleateBoiling.H"
#include "addToRunTimeSelectionTable.H"
#include "myStringOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferModels
{
    defineTypeNameAndDebug(nucleateBoiling, 0);
    addToRunTimeSelectionTable
    (
        heatTransferModel, 
        nucleateBoiling, 
        FSHeatTransferModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatTransferModels::nucleateBoiling::
nucleateBoiling
(
    const objectRegistry& objReg,
    const dictionary& dict,
    const FSPair& FSPair,
    const wordList& regions
)
:
    heatTransferModel
    (
        objReg,
        dict,
        FSPair,
        regions
    ),
    Twall_(FSPair.structure().Twall()),
    Tf_(FSPair.fluidRef().thermo().T()),
    Tsat_
    (
        //- The interface is always at saturation when doing simulations
        //  with phase change enabled, so use this for Tsat_
        FSPair.fluidRef().mesh().lookupObject<volScalarField>("T.interface")
    ),
    p_(FSPair.fluidRef().mesh().lookupObject<volScalarField>("p")),
    pCrit_(3.5e7), //- Specific to Sodium
    flowFactor_
    (
        flowFactorModel::New
        (
            *this,
            objReg,
            FSPair
        )
    ),
    suppressionFactor_
    (
        suppressionFactorModel::New
        (
            *this,
            objReg,
            FSPair
        )
    ),
    convectionHeatTransfer_
    (
        heatTransferModel::New
        (
            objReg,
            dict.subDict("convectionModel"),
            FSPair,
            //- Use the regions list obtained from the top level  
            //  heatTransferModel dictionary name via myStringOps
            myStringOps::split<word>(this->dictName(), ':')
        )
    ),
    poolBoilingHeatTransfer_
    (
        poolBoilingModel::New
        (
            *this,
            objReg,
            FSPair
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::heatTransferModels::nucleateBoiling::correctHtc
(
    volScalarField& htc
) const
{    
    //- Correct htc by first computing its convective component. Now htc
    //  consists only of the purely convective component
    convectionHeatTransfer_->correctHtc(htc);

    //- Correct the convection enhancement flow factor
    flowFactor_->correct();
    const scalarField& F(flowFactor_->flowFactor());

    //- Correct the nucleate boiling suppression factor
    suppressionFactor_->correct();
    const scalarField& S(suppressionFactor_->suppressionFactor());

    //- Correct the pool boiling heat transfer coefficient
    poolBoilingHeatTransfer_->correct();
    const scalarField& htcPB(poolBoilingHeatTransfer_->htc());

    forAll(cellList_, i)
    {
        const label& celli(cellList_[i]);

        //- Recall that htc[celli] on the right hand side consists only of
        //  the convective component
        htc[celli] = F[i]*htc[celli]+S[i]*htcPB[i];
    }
}


// ************************************************************************* //
