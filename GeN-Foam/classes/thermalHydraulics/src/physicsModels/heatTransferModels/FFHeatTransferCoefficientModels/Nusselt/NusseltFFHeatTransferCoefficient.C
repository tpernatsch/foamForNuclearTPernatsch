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

#include "FFPair.H"
#include "NusseltFFHeatTransferCoefficient.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace FFHeatTransferCoefficientModels
{
    defineTypeNameAndDebug(Nusselt, 0);
    addToRunTimeSelectionTable
    (
        FFHeatTransferCoefficientModel, 
        Nusselt, 
        FFHeatTransferCoefficientModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FFHeatTransferCoefficientModels::Nusselt::Nusselt
(
    const FFPair& pair,
    const dictionary& dict,
    const objectRegistry& objReg
)
:
    FFHeatTransferCoefficientModel
    (
        pair,
        dict,
        objReg
    ),
    Re_(pair.Re()),
    kappa_(bulkFluid_.kappa()),
    Pr_(bulkFluid_.Pr()),//(bulkFluid_.Pr()),//(pair.PrContinuous()),
    Dh_(bulkFluid_.Dh()),//(bulkFluid_.Dh()),//(pair.DhDispersed()),
    A_(dict.get<scalar>("const")),
    B_(dict.get<scalar>("coeff")),
    C_(dict.get<scalar>("expRe")),
    D_(dict.get<scalar>("expPr")),
    usePeclet_(C_ == D_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::FFHeatTransferCoefficientModels::Nusselt::value
(
    const label& celli
) const
{
    //- I am creating a scalar on return to (hopefully) force Return Value
    //  Optimizations (RVOs, C++ performance stuff)
    scalar Dhi(max(Dh_[celli], 1e-4));
    if (B_ != 0)
    {
        if (usePeclet_)
            return
                scalar
                ( 
                    (kappa_[celli]/Dhi)*
                    (A_ + B_*pow(Re_[celli]*Pr_[celli], C_))
                );
        else
            return 
                scalar
                (
                    (kappa_[celli]/Dhi)*
                    (A_ + B_*pow(Re_[celli], C_)*pow(Pr_[celli], D_))
                );
    }
    else
        return scalar((kappa_[celli]/Dhi)*A_);

}


// ************************************************************************* //
