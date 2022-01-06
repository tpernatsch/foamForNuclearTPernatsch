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
    Foam::structures::fixedTemperature

Description
    Describes a structure with a constant surface temperature

Authors
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)
    Carlo Fiorina (carlo.fiorina@epfl.ch / carlo.fiorina@outlook.com)

SourceFiles
    fixedTemperature.C

\*---------------------------------------------------------------------------*/

#ifndef fixedTemperature_H
#define fixedTemperature_H

#include "powerModel.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace powerModels
{

/*---------------------------------------------------------------------------*\
                           Class fixedTemperature Declaration
\*---------------------------------------------------------------------------*/

class fixedTemperature
:
    public powerModel
{
protected :

    //- Surface temperature
    mutable volScalarField T_;

    //- Optional table to describe any possible dependence of the temperature
    //  on time
    PtrList<Function1<scalar>> timeProfile_;

    //- Optional offset time
    List<scalar> t0_;

    //- True if a timeProfilePtr_ exists
    List<bool> timeDependent_;

public:

    //- Runtime type information
    TypeName("fixedTemperature");

    //- Constructors

        fixedTemperature
        (
            const structure& structure,
            const dictionary& dicts
        );

        //- Destructor
        virtual ~fixedTemperature();

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