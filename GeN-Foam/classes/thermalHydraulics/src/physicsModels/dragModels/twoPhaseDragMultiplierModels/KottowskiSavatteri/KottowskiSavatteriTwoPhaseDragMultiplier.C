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

#include "KottowskiSavatteriTwoPhaseDragMultiplier.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace twoPhaseDragMultiplierModels
{
    defineTypeNameAndDebug(KottowskiSavatteri, 0);
    addToRunTimeSelectionTable
    (
        twoPhaseDragMultiplierModel, 
        KottowskiSavatteri, 
        twoPhaseDragMultiplierModels
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseDragMultiplierModels::KottowskiSavatteri::KottowskiSavatteri
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
    )
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::twoPhaseDragMultiplierModels::KottowskiSavatteri::phi2
(
    const label& celli
) const
{
    if (onePhase(celli))
    {
        return 1.0;
    }
    //- The correlation is valid only for 7e-2 < X < 30
    scalar log10X(log10(min(max(mFluidPtr_->XLM()[celli], 0.07), 30))); 
    return 
        pow
        (
            10,
            2.0*
            (
                0.1046*sqr(log10X)
            -   0.5098*log10X
            +   0.6252
            )
        );
}

// ************************************************************************* //
