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
#include "FFPair.H"
#include "SahaZuberSubCooledBoilingFraction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace subCooledBoilingFractionModels
{
    defineTypeNameAndDebug(SahaZuber, 0);
    addToRunTimeSelectionTable
    (
        subCooledBoilingFractionModel,
        SahaZuber,
        subCooledBoilingFractionModels
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::subCooledBoilingFractionModels::SahaZuber::SahaZuber
(
    const FSPair& pair,
    const dictionary& dict,
    const objectRegistry& objReg
)
:
    subCooledBoilingFractionModel
    (
        pair,
        dict,
        objReg
    )
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::subCooledBoilingFractionModels::SahaZuber::value
(
    const label& celli,
    const scalar& qi
) const
{
    this->setPtr();

    const scalar& Tsati(FFPairPtr_->iT()[celli]);

    //- Minimum liquid tempreature for bubble detachment
    scalar Tldi =
        max
        (
            Tsati
        -   qi*pair_.structureRef().Dh()[celli]/
            (
                0.0065*max
                (
                    pair_.Re()[celli]*pair_.fluidRef().Pr()[celli], 
                    7e4
                )
            ),
            0.0
        ); 
    return 
        min
        (
            max
            (
                (pair_.fluidRef().thermo().T()[celli] - Tldi)/
                max(Tsati-Tldi, 1e-3),
                0.0
            ),
            1
        );
}

// ************************************************************************* //
