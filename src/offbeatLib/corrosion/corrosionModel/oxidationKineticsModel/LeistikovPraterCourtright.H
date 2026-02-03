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
    Foam::LeistikovPraterCourtright

Description
    High temperature oxidation of Zircaloy-4 combining the intermediate
    temperature model of Leistikov with high temperature model of Prater 
    and Courtright as recommended by Schanz.

    G. Schanz, Recommendations and Supporting Information on the Choice of
    Zirconium Oxidation Models in Severe Accident Codes, Technical report
    FZKA 6827, Forschungszentrum Karlsruhe, 2003.

SourceFiles
    LeistikovPraterCourtright.C

\*---------------------------------------------------------------------------*/

#ifndef LeistikovPraterCourtright_H
#define LeistikovPraterCourtright_H

#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class LeistikovPraterCourtright Declaration
\*---------------------------------------------------------------------------*/

class LeistikovPraterCourtright
{
public:

    //- Runtime type information
    ClassName("Leistikov-Prater-Courtright");

    //- Construct from dictionary
    LeistikovPraterCourtright(const dictionary& dict) {}
    
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
