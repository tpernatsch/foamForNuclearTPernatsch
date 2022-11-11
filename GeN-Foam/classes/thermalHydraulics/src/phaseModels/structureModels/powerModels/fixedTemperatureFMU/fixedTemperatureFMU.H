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
    Foam::structures::fixedTemperatureFMU

Description
    Describes a structure with a temperature taken from an FMU

Dictionary entries example (mandatory)
\verbatim
    type            NusseltWallAndHfromFMU;
    HNameFromFMU    HCoupled;
    const   0.0;
    coeff   0.023;
    expRe   0.8;
    expPr   0.4;
    addH    30000;
\endverbatim

Authors
    Carlo Fiorina (carlo.fiorina@epfl.ch / carlo.fiorina@outlook.com)

SourceFiles
    fixedTemperatureFMU.C

\*---------------------------------------------------------------------------*/

#ifndef fixedTemperatureFMU_H
#define fixedTemperatureFMU_H

#include "powerModel.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace powerModels
{

/*---------------------------------------------------------------------------*\
                           Class fixedTemperatureFMU Declaration
\*---------------------------------------------------------------------------*/

class fixedTemperatureFMU
:
    public powerModel
{
protected :

    //- Surface temperature
    mutable volScalarField T_;

    //- Name of temperature from FMU
    word temperatureNameFromFMU_;

public:

    //- Runtime type information
    TypeName("fixedTemperatureFMU");

    //- Constructors

        fixedTemperatureFMU
        (
            structure& structure,
            const dictionary& dicts
        );

        //- Destructor
        virtual ~fixedTemperatureFMU();

    //- Member Functions

        //- In two-phase solvers :
        //  HT = H1*frac1*T1 + H2*frac2*T2
        //  H =  H1*frac1 + H2*frac2
        //- In a mono-phase solver :
        //  HT = H*T
        //  H = H
        virtual void correct
        (
            const volScalarField& HT,
            const volScalarField& H
        ){}

        //- Set provided field with structure surface temperature
        virtual void correctT(volScalarField& T) const;

        //- Turn off power
        virtual void powerOff();

        //- Temperature update based on given profile
        void temperatureUpdate() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace powerModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //