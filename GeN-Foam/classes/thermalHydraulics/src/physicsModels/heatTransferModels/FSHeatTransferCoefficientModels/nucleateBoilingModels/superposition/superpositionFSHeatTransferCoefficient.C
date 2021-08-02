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
#include "superpositionFSHeatTransferCoefficient.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace FSHeatTransferCoefficientModels
{
    defineTypeNameAndDebug(superposition, 0);
    addToRunTimeSelectionTable
    (
        FSHeatTransferCoefficientModel, 
        superposition, 
        FSHeatTransferCoefficientModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FSHeatTransferCoefficientModels::superposition::superposition
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
    htcFCPtr_
    (
        FSHeatTransferCoefficientModel::New
        (
            pair,
            this->subDict("forcedConvection"),
            pair.mesh()
        )
    ),
    htcPBPtr_
    (
        FSHeatTransferCoefficientModel::New
        (
            pair,
            this->subDict("poolBoiling"),
            pair.mesh()
        )
    ),
    FPtr_
    (
        flowEnhancementFactorModel::New
        (
            pair,
            this->subDict("flowEnhancementFactor"),
            pair.mesh()
        )
    ),
    SPtr_
    (
        suppressionFactorModel::New
        (
            pair,
            this->subDict("suppressionFactor"),
            pair.mesh()
        )
    )
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::FSHeatTransferCoefficientModels::superposition::value
(
    const label& celli
) const
{
    return 
        htcFCPtr_->value(celli)*FPtr_->value(celli)
    +   htcPBPtr_->value(celli)*SPtr_->value(celli);
}

// ************************************************************************* //
