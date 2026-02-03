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

#include "EPRI_SLI.H"
#include "scalar.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(EPRI_SLI, 0);

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

EPRI_SLI::EPRI_SLI
(
    const dictionary& dict
)
:
    Li_(readScalar(dict.lookup("Li"))),
    Sn_(readScalar(dict.lookup("Sn"))),
    Fe_(readScalar(dict.lookup("Fe")))
{}


void EPRI_SLI::correctOxideThickness
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
    const scalar dt = dt_seconds / (3600*24);   // Time in days (d)
    const scalar R = 1.987;                     // Ideal gas constant (cal/mol.K)
    
    //- TODO: Include correct hydrogen concentration
    const scalar H2 = 0;                        // Assume zero hydrogen (ppm)
    //- TODO: Include correct heat flux
    const scalar kox = 0.94;
    const scalar QA = kox / max(S0, SMALL) * (Ti - Tb); // Approximate heat flux from temperature difference
    //const scalar QA = 0;
    
    const scalar S_trans = 2e-6;                // Transition thickness (m)
    const scalar C1 = 1e-18*5.867e10;           // (m3/day)
    const scalar C2 = 1e-6*7.619e6;             // (m/day)
    const scalar C_Li = 0.65;
    const scalar C_Fe = 2e-2;
    const scalar C_phi = 1.2e-4;
    const scalar P0 = 0.24;
    const scalar Q1R = 33662.7 / R;
    const scalar Q2R = max(24825 + (9135.6 - 24825)*H2, 9135.6) / R;
    
    const scalar F_Li1 = exp(C_Li*(0.12*Li_ - 23.0*Li_ / Tb));
    const scalar F_Li2 = exp(C_Li*(0.17*Li_ - 20.4*Li_ / Tb));
    const scalar F_Fe = C_Fe*Fe_;
    const scalar F_Sn = min(1.25*(Sn_-1.38) + 1, 0.75*(Sn_-1.38) + 1);
    const scalar F_phi = C_phi*pow(phi, P0);
    const scalar F_QA = 1 + 8.81e-8*QA;
    const scalar F_H2 = (H2 > 400 ? 0.699*log(H2/400) : 0);
    
    const scalar k1 = C1*F_Li1*(1 + F_Fe);
    const scalar k2 = C2*F_Li2*F_Sn*F_QA*(1 + F_H2 + F_Fe + F_phi);

    // Compute oxide thickness
    if (S0 < S_trans)
    {
        scalar dt_trans = (pow(S_trans, 3.0) - pow(S0, 3.0))/(C1*exp(-Q1R/Ti));
    
        if (dt_trans >= dt)
        {
            // Pre-transition regime
            S = pow(k1*exp(-Q1R/Ti)*dt + pow(S0, 3.0), 1.0/3.0);
        }
        else
        {
            // Transitioning
            scalar dt_post = dt - dt_trans;
            
            S = S_trans + k2*exp(-Q2R/Ti)*dt_post;            
        }
    }
    else
    {
        // Post-transition regime
        S = S0 + k2*exp(-Q2R/Ti)*dt;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
