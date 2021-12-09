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
    Foam::structures::timeDependentPower

Description
    Describes a generic structure with a time dependent internal power production.
    Solves an energy equation for the structure, updating the average
    structure temperature (assumed equal to the structure surface temperature
    for fluid heat transfer purposes).

Author
    Carlo Fiorina (carlo.fiorina@epfl.ch / carlo.fiorina@outlook.com)

SourceFiles
    timeDependentPower.C

\*---------------------------------------------------------------------------*/

#ifndef timeDependentPower_H
#define timeDependentPower_H

#include "powerModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Forward declaration
class structure;

namespace powerModels
{

/*---------------------------------------------------------------------------*\
                           Class timeDependentPower Declaration
\*---------------------------------------------------------------------------*/

class timeDependentPower
:
    public powerModel
{
protected :

    //- Volumetric power
    volScalarField powerDensity_;

    //- Time dependent power
    autoPtr<Function1<scalar>> timeDependentPowerPtr_;    

    //- Structure temperature
    volScalarField T_;

    //- Volumetric heat capacity of the structure
    volScalarField alphaRhoCp_;

public:

    //- Runtime type information
    TypeName("timeDependentPower");

    //- Constructors

        timeDependentPower
        (
            const structure& structure,
            const dictionary& dicts
        );

        //- Destructor
        virtual ~timeDependentPower();

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
        );

        //- Set provided field with structure surface temperature
        virtual void correctT(volScalarField& T) const;

        //- Power off
        virtual void powerOff();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace powerModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //