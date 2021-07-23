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

Class

Description
    Calculate the temperature of onset of nucleate boiling according to the
    model by N. Basu, G.R. Warrier, and V.K. Dhir, "Onset of Nucleate Boiling 
    and Active Nucleation Site Density During Subcooled Flow Boiling," 
    J. Heat Transfer, 124, 717-728, 2002.

Author
    Stefan Radman (sradman@protonmail.com / stefan.radman@epfl.ch)

SourceFiles
    BasuTONB.C

\*---------------------------------------------------------------------------*/

#ifndef BasuTONB_H
#define BasuTONB_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "TONBModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace TONBModels
{

/*---------------------------------------------------------------------------*\
            Class Basu Declaration
\*---------------------------------------------------------------------------*/

class Basu
:
    public TONBModel
{
protected:

    //- Liquid surface tension
    scalar sigma_;

    //- Liquid-structure interface contact angle
    scalar contactAngle_;

    //- Pre-computed constant
    scalar A_;

public:

    //- Runtime type information
    TypeName("Basu");

    // Constructors

        Basu
        (
            const FSPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~Basu(){}

    // Member Functions

        //- Return TONB
        virtual scalar value
        (
            const label& celli, 
            const scalar& htc2pFCi
        ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace TONBModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
