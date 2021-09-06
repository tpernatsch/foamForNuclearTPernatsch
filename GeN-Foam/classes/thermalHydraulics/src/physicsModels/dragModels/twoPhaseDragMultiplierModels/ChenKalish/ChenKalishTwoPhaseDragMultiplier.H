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
    Foam::ChenKalish

Description
    Model by Chen and Kalish for the multiplier, specifically for
    potassium two-phase flows:

    phi2 = e^(2*( 0.0867*(logX)^2 - 0.518*logX + 1.59 ))

    log is the natural logarithm, while X2 is the well-known 
    Lockhart-Martinelli parameter and is defined as:

    X2 = (mu_I/mu_O)^0.2 * (flowQuality_I/flowQuality_O)^1.8 * rho_O/rho_I

    Chen, J.C., Kalish, S., 1970. "An experimental investigation of two-phase
    pressure drop for potassium with and without net vaporization." In: 4th 
    Int. Heat transfer Conf, Vol. 6.

Author
    Stefan Radman (sradman@protonmail.com / stefan.radman@epfl.ch)

SourceFiles
    ChenKalish.C

\*---------------------------------------------------------------------------*/

#ifndef ChenKalishTwoPhaseDragMultiplier_H
#define ChenKalishTwoPhaseDragMultiplier_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "twoPhaseDragMultiplierModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace twoPhaseDragMultiplierModels
{

/*---------------------------------------------------------------------------*\
                    Class ChenKalish Declaration
\*---------------------------------------------------------------------------*/

class ChenKalish
:
    public twoPhaseDragMultiplierModel
{
public:

    //- Runtime type information
    TypeName("ChenKalish");

    // Constructors

        ChenKalish
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~ChenKalish(){}

    // Member Functions

        //- Return multiplier
        virtual scalar phi2(const label& celli) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace twoPhaseDragMultiplierModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
