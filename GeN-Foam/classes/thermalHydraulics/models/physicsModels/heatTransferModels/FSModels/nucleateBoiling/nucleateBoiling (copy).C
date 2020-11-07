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
    convectionHeatTransferModel_
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
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::heatTransferModels::nucleateBoiling::correctHtc
(
    volScalarField& htc
) const
{    
    scalar C(6.9);
    scalar m(0.12);
    scalar n(0.7);

    //- Correct htc by first computing its convective component
    convectionHeatTransferModel_->correctHtc(htc);

    //- Correct the convection enhancement flow factor
    flowFactor_->correct();
    const scalarField& F(flowFactor_->flowFactor());

    forAll(cellList_, i)
    {
        const label& celli(cellList_[i]);
        
        //- Forced convection component. Recall that htc was set to the 
        //  convective component only by the correctHtc called by
        //  convectionHeatTransferModel
        scalar htcFC(htc[celli]);
        
        //- Reduced pressure
        scalar pR(p_[celli]/pCrit_);

        scalar deltaT(max(Twall_[celli]-Tf_[celli], 1e-9));

        //- Pool boiling HTC from Shah (specifically for sodium)
        scalar htcPB
        (
            pow
            (
                C*pow(pR, m)*pow(deltaT, n),
                1.0/(1.0-n)
            )
        );

        //- Suppression factor
        scalar S(max(Tf_[celli]-Tsat_[celli], 0)/deltaT);

        //- Forced convection enhancement factor
        htc[celli] = F[i]*htcFC+S*htcPB;
    }

    /*
    if (B_ != 0)
    {
        if (usePeclet_)
        {
            forAll(cellList_, i)
            {
                label celli(cellList_[i]);
                htc[celli] = 
                    (kappa[celli]/Dh_[celli])*
                    (A_ + B_*pow(Re_[celli]*Pr_[celli], C_));
            }
        }
        else
        {
            forAll(cellList_, i)
            {
                label celli(cellList_[i]);
                htc[celli] = 
                    (kappa[celli]/Dh_[celli])*
                    (A_ + B_*pow(Re_[celli], C_)*pow(Pr_[celli], D_));
            }
        }
    }
    else
    {
        forAll(cellList_, i)
        {
            label celli(cellList_[i]);
            htc[celli] = (kappa[celli]/Dh_[celli])*A_;
        }
    }
    */
}


// ************************************************************************* //
