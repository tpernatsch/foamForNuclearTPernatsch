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

#include "MATPRO_CORROS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(MATPRO_CORROS_PWR, 0);

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void MATPRO_CORROS_PWR::correctOxideThickness
(
    scalar& S,
    const scalar S0,
    const scalar Ti,
    const scalar Tb,
    const scalar phi,
    const scalar dt_second
) const
{
    // Correlation parameters
    const scalar A = 1.203e2*exp(-7.118e-3*Tb);     // PWR conditions
    const scalar C1 = 4.976e-9*A;                   // [m3/day]
    const scalar C2 = 82.88*A;	                    // [m/day]
    const scalar Q1R = 15660;                       // [1/K]
    const scalar Q2R = 14080;                       // [1/K]

    const scalar dt = dt_second / 86400;            // Time in days

    // Transition thickness
    scalar S_trans = 7.749e-6*exp(-790.0/Ti);

    // Compute oxide thickness
    if( S0 < S_trans )
    {
        // Get time of transition
        scalar dt_trans = (pow(S_trans, 3) - pow(S0, 3)) / (C1*exp(-Q1R/Ti));
        
        if (dt_trans > dt)
        {
            // Pre-transition regime
            S = pow(C1*exp(-Q1R/Ti)*dt + pow(S0, 3.0), 1.0/3.0);
        }
        else
        {
            // Transition
            scalar dt_post = dt - dt_trans;
            
            S = S_trans + C2*exp(-Q2R/Ti)*dt_post;
        }
    }
    else
    {
        // Post transition regime
        S = S0 + C2*exp(-Q2R/Ti)*dt;
    }
}


void MATPRO_CORROS_BWR::correctOxideThickness
(
    scalar& S,
    const scalar S0,
    const scalar Ti,
    const scalar Tb,
    const scalar phi,
    const scalar dt_second
) const
{
    // Correlation parameters
    const scalar A = 4.84e5*exp(-1.945e-2*Tb);      // BWR conditions
    const scalar C1 = 4.976e-9*A;                   // [m3/day]
    const scalar C2 = 8.288e1*A;	                // [m/day]
    const scalar Q1R = 15660;                       // [1/K]
    const scalar Q2R = 14080;                       // [1/K]

    const scalar dt = dt_second / 86400;            // Time in days

    // Transition thickness
    scalar S_trans = 7.749e-6*exp(-790.0/Ti);

    // Compute oxide thickness
    if( S0 < S_trans )
    {
        // Get time of transition
        scalar dt_trans = (pow(S_trans, 3) - pow(S0, 3)) / (C1*exp(-Q1R/Ti));
        
        if (dt_trans > dt)
        {
            // Pre-transition regime
            S = pow(C1*exp(-Q1R/Ti)*dt + pow(S0, 3.0), 1.0/3.0);
        }
        else
        {
            // Transition
            scalar dt_post = dt - dt_trans;
            
            S = S_trans + C2*exp(-Q2R/Ti)*dt_post;
        }
    }
    else
    {
        // Post transition regime
        S = S0 + C2*exp(-Q2R/Ti)*dt;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
