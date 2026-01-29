/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "LeistikovPraterCourtright.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(LeistikovPraterCourtright, 0);

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void LeistikovPraterCourtright::correctOxideThickness
(
    scalar& S,
    const scalar S0,
    const scalar Ti,
    const scalar Tb,
    const scalar phi,
    const scalar dt
) const
{
    // Correlation parameters
    const scalar A1  = 7.82e-6;
    const scalar Q1R = 20214;
    const scalar T1 = 1800;
    const scalar A2  = 2.98116e-3;
    const scalar Q2R = 28420;
    const scalar T2 = 1900;

    scalar As;
    scalar QsR;

    // Leistikov correlation
    if ( Ti < T1 )
    {
        QsR = Q1R;
        As = A1;
    }
    // Interpolation
    else if (Ti < T2)
    {
        scalar y1 = log(A1) - Q1R/T1;
        scalar y2 = log(A2) - Q2R/T2;

        QsR = (y1 - y2) / (1/T2 - 1/T1);
        As = exp(y1 + QsR/T1);
    }
    // Prater-Courtright correlation
    else
    {
        QsR = Q2R;
        As = A2;
    }

    S = sqrt( As*exp(-QsR/Ti)*dt + pow(S0, 2.0) );
}

} // End namespace Foam

// ************************************************************************* //
