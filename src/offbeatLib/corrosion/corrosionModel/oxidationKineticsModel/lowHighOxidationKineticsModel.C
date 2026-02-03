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

#include "lowHighOxidationKineticsModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Lo, class Hi>
lowHighOxidationKineticsModel<Lo, Hi>::lowHighOxidationKineticsModel
(
    const dictionary& dict
)
:
    oxidationKineticsModel(dict),
    lowTModel_(dict),
    highTModel_(dict)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Lo, class Hi>
lowHighOxidationKineticsModel<Lo, Hi>::~lowHighOxidationKineticsModel()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


template<class Lo, class Hi>
void lowHighOxidationKineticsModel<Lo, Hi>::correctOxideThickness
(
    scalarField& S,
    const scalarField& S0,
    const scalarField& Ti,
    const scalarField& Tb,
    const scalarField& phi,
    const scalar dt
) const
{
    forAll(S, i)
    {
        // Normal operating conditions
        lowTModel_.correctOxideThickness
        (
            S[i], S0[i], Ti[i], Tb[i], phi[i], dt
        );
        
        if (Ti[i] < lowTModel_.upperLimit())
        {
            // Normal operating conditions
            lowTModel_.correctOxideThickness
            (
                S[i], S0[i], Ti[i], Tb[i], phi[i], dt
            );
        }
        else
        {
            scalar sHi=S[i];
            scalar sLo = S[i];
            
            // Normal operating conditions
            lowTModel_.correctOxideThickness
            (
                sLo, S0[i], lowTModel_.upperLimit(), min(Tb[i], lowTModel_.upperLimit()), phi[i], dt
            );
        
            // High temperature model:
            // Normal operating conditions
            highTModel_.correctOxideThickness
            (
                sHi, S0[i], Ti[i], Tb[i], phi[i], dt
            );
            
            S[i] = max(sLo, sHi);
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
