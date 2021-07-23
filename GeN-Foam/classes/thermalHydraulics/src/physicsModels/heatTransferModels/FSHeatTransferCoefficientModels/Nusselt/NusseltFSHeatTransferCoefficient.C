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
#include "NusseltFSHeatTransferCoefficient.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace FSHeatTransferCoefficientModels
{
    defineTypeNameAndDebug(Nusselt, 0);
    addToRunTimeSelectionTable
    (
        FSHeatTransferCoefficientModel, 
        Nusselt, 
        FSHeatTransferCoefficientModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FSHeatTransferCoefficientModels::Nusselt::Nusselt
(
    const FSPair& pair,
    const dictionary& dict,
    const objectRegistry& objReg
)
:
    FSHeatTransferCoefficientModel
    (
        pair,
        dict,
        objReg
    ),
    Re_(pair.Re()),
    kappa_(pair.fluidRef().kappa()),
    Pr_(pair.fluidRef().Pr()),
    Dh_(pair.fluidRef().Dh()),
    A_(dict.get<scalar>("const")),
    B_(dict.get<scalar>("coeff")),
    C_(dict.get<scalar>("expRe")),
    D_(dict.get<scalar>("expPr")),
    usePeclet_(C_ == D_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::FSHeatTransferCoefficientModels::Nusselt::value
(
    const label& celli
) const
{
    //- I am creating a scalar on return to (hopefully) force Return Value
    //  Optimizations (RVOs, C++ performance stuff)
    if (B_ != 0)
    {
        if (usePeclet_)
            return
                scalar
                ( 
                    (kappa_[celli]/Dh_[celli])*
                    (A_ + B_*pow(Re_[celli]*Pr_[celli], C_))
                );
        else
            return 
                scalar
                (
                    (kappa_[celli]/Dh_[celli])*
                    (A_ + B_*pow(Re_[celli], C_)*pow(Pr_[celli], D_))
                );
    }
    else
        return scalar((kappa_[celli]/Dh_[celli])*A_);
}


// ************************************************************************* //
