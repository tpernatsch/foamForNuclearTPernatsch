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
    Foam::EPRI_KWU_CE

Description
    EPRI/KWU/C-E low-temperature oxidation model, taken from EPRI NP-2789.
    
    Waterside Corrosion of Zircaloy Fuel Rods, EPRI NP-2789, December 1982.

SourceFiles
    EPRI_KWU_CE.C

\*---------------------------------------------------------------------------*/

#ifndef EPRI_KWU_CE_H
#define EPRI_KWU_CE_H

#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class EPRI_KWU_CE Declaration
\*---------------------------------------------------------------------------*/

class EPRI_KWU_CE
{
public:

    //- Runtime type information
    ClassName("EPRI-KWU-CE");

    //- Construct from dictionary
    EPRI_KWU_CE(const dictionary& dict) {}
    
    // Member Functions

        //- Correlation lower temperature limit
        inline scalar lowerLimit() const
        {
            return 500;
        }

        //- Correlation upper temperature limit
        inline scalar upperLimit() const
        {
            return 673;
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
