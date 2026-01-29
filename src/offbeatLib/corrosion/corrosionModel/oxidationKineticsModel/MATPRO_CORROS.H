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
    Foam::MATPRO_CORROS

Description
    MATPRO 11 low-temperature oxidation model based on CORROS subcode.

    MATPRO — VERSION 11
    A HANDBOOK OF MATERIALS PROPERTIES FOR USE IN THE ANALYSIS OF LIGHT
    WATER REACTOR FUEL ROD BEHAVIOR, NUREG/CR-0497, February 1979.

SourceFiles
    MATPRO_CORROS.C

\*---------------------------------------------------------------------------*/

#ifndef MATPRO_CORROS_H
#define MATPRO_CORROS_H

#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class MATPRO_CORROS Declaration
\*---------------------------------------------------------------------------*/

class MATPRO_CORROS_PWR
{
public:

    //- Runtime type information
    ClassName("MATPRO-CORROS-PWR");

    //- Construct from dictionary
    MATPRO_CORROS_PWR(const dictionary& dict) {}
    
    // Member Functions

        //- Correlation lower temperature limit
        inline scalar lowerLimit() const
        {
            return 523.15;
        }

        //- Correlation upper temperature limit
        inline scalar upperLimit() const
        {
            return 673.15;
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

class MATPRO_CORROS_BWR
{
public:

    //- Runtime type information
    ClassName("MATPRO-CORROS-BWR");

    //- Construct from dictionary
    MATPRO_CORROS_BWR(const dictionary& dict) {}
    
    // Member Functions

        //- Correlation lower temperature limit
        inline scalar lowerLimit() const
        {
            return 523.15;
        }

        //- Correlation upper temperature limit
        inline scalar upperLimit() const
        {
            return 673.15;
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
