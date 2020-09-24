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

#include "NusseltReynoldsPrandtlPower.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferModels
{
    defineTypeNameAndDebug(NusseltReynoldsPrandtlPower, 0);
    addToRunTimeSelectionTable
    (
        heatTransferModel, 
        NusseltReynoldsPrandtlPower, 
        FFHeatTransferModels
    );
    addToRunTimeSelectionTable
    (
        heatTransferModel, 
        NusseltReynoldsPrandtlPower, 
        FSHeatTransferModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatTransferModels::NusseltReynoldsPrandtlPower::
NusseltReynoldsPrandtlPower
(
    const objectRegistry& objReg,
    const dictionary& dict,
    const FFPair& FFPair,
    const word nameBulk,
    const word nameInterface
)
:
    heatTransferModel
    (
        objReg,
        dict,
        FFPair,
        nameBulk,
        nameInterface
    ),
    bulkFluid_
    (
        (
            FFPair.fluid1().name() == nameBulk
        ) ?
        FFPair.fluid1() :
        FFPair.fluid2()
    ),
    Re_(FFPair.Re()),
    Pr_(FFPair.PrContinuous()),
    Dh_(FFPair.DhDispersed()),
    A_(dict.get<scalar>("const")),
    B_(dict.get<scalar>("coeff")),
    C_(dict.get<scalar>("expRe")),
    D_(dict.get<scalar>("expPr")),
    usePeclet_(C_ == D_)
{
    if 
    (
            FFPair.fluid1().name() != nameBulk 
        and FFPair.fluid2().name() != nameBulk
    )
    {
        FatalErrorInFunction
            << "Phase " << nameBulk << " not found!"
            << exit(FatalError);
    }
}


Foam::heatTransferModels::NusseltReynoldsPrandtlPower::
NusseltReynoldsPrandtlPower
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
    bulkFluid_(FSPair.fluidRef()),
    Re_(FSPair.Re()),
    Pr_(FSPair.fluidRef().Pr()),
    Dh_(FSPair.structure().Dh()),
    A_(dict.get<scalar>("const")),
    B_(dict.get<scalar>("coeff")),
    C_(dict.get<scalar>("expRe")),
    D_(dict.get<scalar>("expPr")),
    usePeclet_(C_ == D_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::heatTransferModels::NusseltReynoldsPrandtlPower::correctHtc
(
    volScalarField& htc
) const
{    
    //- NEVER use a function that returns a tmp to init a const ref, do this
    //  instead
    tmp<volScalarField> tkappa(bulkFluid_.thermo().kappa());
    volScalarField& kappa = tkappa.ref();

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
}


// ************************************************************************* //
