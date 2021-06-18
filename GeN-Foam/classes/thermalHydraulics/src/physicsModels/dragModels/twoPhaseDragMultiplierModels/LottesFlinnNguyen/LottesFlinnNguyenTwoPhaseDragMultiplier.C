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

#include "LottesFlinnNguyenTwoPhaseDragMultiplier.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace twoPhaseDragMultiplierModels
{
    defineTypeNameAndDebug(LottesFlinnNguyen, 0);
    addToRunTimeSelectionTable
    (
        twoPhaseDragMultiplierModel, 
        LottesFlinnNguyen, 
        twoPhaseDragMultiplierModels
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseDragMultiplierModels::LottesFlinnNguyen::LottesFlinnNguyen
(
    const fvMesh& mesh,
    const dictionary& dict,
    const objectRegistry& objReg
)
:
    twoPhaseDragMultiplierModel
    (
        mesh,
        dict,
        objReg
    ),
    exp_(dict.lookupOrDefault<scalar>("exp", 2.0))
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::twoPhaseDragMultiplierModels::LottesFlinnNguyen::phi2
(
    const label& celli
) const
{
    //- If it is "almost" single-phase flow, return 1.0
    return 
        (mFluidPtr_->normalized()[celli] > 0.999) ?
        1.0 :
        pow(1.0-pow(1.0+pow(mFluidPtr_->XLM()[celli], 0.8), -0.378), -exp_);
}

// ************************************************************************* //
