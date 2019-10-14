/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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
    Foam::saturationModels::constantTemperature

Description
    Constant saturation pressure and temperature.

SourceFiles
    constantTemperature.C

\*---------------------------------------------------------------------------*/

#ifndef constantTemperature_H
#define constantTemperature_H

#include "saturationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace saturationModels
{

/*---------------------------------------------------------------------------*\
                         Class constantTemperature Declaration
\*---------------------------------------------------------------------------*/

class constantTemperature
:
    public saturationModel
{
protected:

    //- Constant saturation temperature
    dimensionedScalar Tsat_;

public:

    TypeName("constantTemperature");

    //- Constructors

        //- Construct from a dictionary
        constantTemperature
        (
            const dictionary& dict, 
            const fluid& fluid1,
            const fluid& fluid2
        );

        //- Destructor
        virtual ~constantTemperature(){}

    //- Member Functions

        //- Saturation temperature at given pressure
        virtual tmp<volScalarField> Tsat
        (
            const volScalarField& p
        ) const;

        //- Saturation pressure at given temperature
        virtual tmp<volScalarField> pSat
        (
            const volScalarField& T
        ) const;

        //- Natural log of the saturation pressure at given temperature
        virtual tmp<volScalarField> lnPSat
        (
            const volScalarField& T
        ) const;   

        //- dPressureSat/dTemperature at given temperature
        virtual tmp<volScalarField> pSatPrime
        (
            const volScalarField& T
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace saturationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
