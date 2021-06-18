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

#include "Kaiser88TwoPhaseDragMultiplier.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace twoPhaseDragMultiplierModels
{
    defineTypeNameAndDebug(Kaiser88, 0);
    addToRunTimeSelectionTable
    (
        twoPhaseDragMultiplierModel, 
        Kaiser88, 
        twoPhaseDragMultiplierModels
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseDragMultiplierModels::Kaiser88::Kaiser88
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

Foam::scalar Foam::twoPhaseDragMultiplierModels::Kaiser88::phi2
(
    const label& celli
) const
{
    if (onePhase(celli))
        return 1.0;
    //- I am limiting this for 7e-2 < X < 30 like Kottowski-Savatteri
    //  out of consistency
    scalar logSqrtX(log(sqrt(min(max(mFluidPtr_->XLM()[celli], 0.07), 30))));
    return 
        exp
        (
            2.0*
            (
                1.48
            -   1.05*logSqrtX
            +   0.09*sqr(logSqrtX)
            )
        );
}

// ************************************************************************* //
