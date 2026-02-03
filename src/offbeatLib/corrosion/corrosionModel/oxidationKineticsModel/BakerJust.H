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

Class
    Foam::BakerJust

Description
    Baker-Just high temperature oxidation model

SourceFiles
    BakerJust.C

\*---------------------------------------------------------------------------*/

#ifndef BakerJust_H
#define BakerJust_H

#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class EPRI_KWU_CE Declaration
\*---------------------------------------------------------------------------*/

class BakerJust
{
public:

    //- Runtime type information
    ClassName("Baker-Just");

    //- Construct from dictionary
    BakerJust(const dictionary& dict) {}
    
    // Member Functions

        //- Correlation lower temperature limit
        inline scalar lowerLimit() const
        {
            return 673;
        }

        //- Correlation upper temperature limit
        inline scalar upperLimit() const
        {
            return GREAT;
        }

        //- Correct the oxide layer thickness [m]
        void correctOxideThickness
        (
            scalar& S,
            const scalar S0,
            const scalar Ti,
            const scalar Tb,
            const scalar phi,
            const scalar dt
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
