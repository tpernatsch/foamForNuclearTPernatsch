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

#include "EPRI_KWU_CE.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(EPRI_KWU_CE, 0);

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void EPRI_KWU_CE::correctOxideThickness
(
    scalar& S,
    const scalar S0,
    const scalar Ti,
    const scalar Tb,
    const scalar phi,
    const scalar dt_seconds
) const
{
    // Correlation parameters
    const scalar dt = dt_seconds / (3600*24);                           // Time in days
    const scalar C1 = 1e-18*6.3e9;										// [m3/day]
    const scalar C2 = 1e-6*(8.04e7 + 2.59e8*pow(7.46e-15*phi,0.24));	// [m/day]
    const scalar Q1R = 32289 / 1.9872;      							// [K], R = 1.9872 cal/mol.K
    const scalar Q2R = 27354 / 1.9872;									// [K], R = 1.9872 cal/mol.K
    const scalar Q3R = 10763 / 1.9872;									// [K], R = 1.9872 cal/mol.K

    const scalar S_trans = 21.4*exp(-Q3R/Ti - 0.0117*Ti);               // Transition thickness
    

    // Compute oxide thickness
    if (S0 < S_trans)
    {
        scalar dt_trans = (pow(S_trans, 3.0) - pow(S0, 3.0))/(C1*exp(-Q1R/Ti));
    
        if (dt_trans >= dt)
        {
            // Pre-transition regime
            S = pow(C1*exp(-Q1R/Ti)*dt + pow(S0, 3.0), 1.0/3.0);
        }
        else
        {
            // Transitioning
            scalar dt_post = dt - dt_trans;
            
            S = S_trans + C2*exp(-Q2R/Ti)*dt_post;            
        }
    }
    else
    {
        // Post-transition regime
        S = S0 + C2*exp(-Q2R/Ti)*dt;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
