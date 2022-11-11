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
    Foam::structures::fixedPower

Description
    Describes a generic structure with a constant or time depdendent internal 
    power production.
    Solves an energy equation for the structure, updating the average
    structure temperature (assumed equal to the structure surface temperature
    for fluid heat transfer purposes).

Author
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)
    Carlo Fiorina (carlo.fiorina@epfl.ch / carlo.fiorina@outlook.com)

SourceFiles
    fixedPower.C

\*---------------------------------------------------------------------------*/

#ifndef fixedPower_H
#define fixedPower_H

#include "powerModel.H"
#include "Function1.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Forward declaration
class structure;

namespace powerModels
{

/*---------------------------------------------------------------------------*\
                           Class fixedPower Declaration
\*---------------------------------------------------------------------------*/

class fixedPower
:
    public powerModel
{
protected :

    //- Volumetric power
    volScalarField powerDensity_;

    //- Optional table to describe any possible dependence of the power
    //  on time
    PtrList<Function1<scalar>> timeProfile_;

    //- Optional offset time
    List<scalar> t0_;

    //- True if a timeProfilePtr_ exists
    List<bool> timeDependent_;

    //- Structure temperature
    volScalarField T_;

    //- Volumetric heat capacity of the structure
    volScalarField alphaRhoCp_;

public:

    //- Runtime type information
    TypeName("fixedPower");

    //- Constructors

        fixedPower
        (
            structure& structure,
            const dictionary& dicts
        );

        //- Destructor
        virtual ~fixedPower();

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

        //- Power update based on given profile
        void powerUpdate();

        //- Access

            //-
            //bool timeDependent() const
            //{
            //    return timeDependent_;
            //}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace powerModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //