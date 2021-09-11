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
#include "BasuTONB.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace TONBModels
{
    defineTypeNameAndDebug(Basu, 0);
    addToRunTimeSelectionTable
    (
        TONBModel,
        Basu,
        TONBModels
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::TONBModels::Basu::Basu
(
    const FSPair& pair,
    const dictionary& dict,
    const objectRegistry& objReg
)
:
    TONBModel
    (
        pair,
        dict,
        objReg
    ),
    sigma_(dict.get<scalar>("surfaceTension")),
    contactAngle_
    (
        dict.get<scalar>("contactAngle")*
        2.0*constant::mathematical::pi/360.0
    ),
    A_
    (
        2.0*sigma_/sqr((1.0-exp(-pow(contactAngle_, 3)-0.5*contactAngle_)))
    )
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::TONBModels::Basu::value
(
    const label& celli,
    const scalar& htc2pFCi
) const
{
    const scalar& Tli(pair_.fluidRef().T()[celli]);
    const scalar& Tsati(Tsat_[celli]);

    if (Tli < Tsati)
    {
        this->setPtrs();
        
        scalar deltaTONBsat
        (
            A_*htc2pFCi*Tsati/
            (
                otherFluidPtr_->rho()[celli]*
                ((*LPtr_)[celli])*
                pair_.fluidRef().kappa()[celli]
            )
        );

        return  
        (
            Tli + 0.25*
            sqr
            (
                sqrt(deltaTONBsat)
            +   sqrt(deltaTONBsat+4.0*(Tsati-Tli))
            )
        );
    }
    else 
        //- if Tli > Tsati then whatever I return is fundamentally irrelevant
        //  due to the over-arching code structure of multiRegimeBoiling
        //  that uses this class
        return 0.0;
}

// ************************************************************************* //
