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

#include "LockhartMartinelli.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace twoPhaseDragMultiplierModels
{
    defineTypeNameAndDebug(LockhartMartinelli, 0);
    addToRunTimeSelectionTable
    (
        twoPhaseDragMultiplierModel, 
        LockhartMartinelli, 
        twoPhaseDragMultiplierModels
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseDragMultiplierModels::LockhartMartinelli::LockhartMartinelli
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
    C_
    (
        dict.lookupOrDefault<scalar>("C", 20)
    )
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::twoPhaseDragMultiplierModels::LockhartMartinelli::phi2
(
    const label& celli
) const
{
    if (onePhase(celli))
        return 1.0;
    const scalar& X(mFluidPtr_->XLM()[celli]);
    return 
        scalar
        (
            min(1.0 + C_/X + 1.0/sqr(X), maxPhi2_)
        );
}

// ************************************************************************* //
