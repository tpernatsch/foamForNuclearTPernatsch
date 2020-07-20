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
    Foam::dragModels::Rehme

Description
    Rheme correlation for pressure drop

    @ARTICLE{Rehme197315,
    author={Rehme, Klaus},
    title={PRESSURE DROP CORRELATIONS FOR FUEL ELEMENT SPACERS.},
    journal={Nuclear Technology},
    year={1973},
    volume={17},
    number={1},
    pages={15-23},
    doi={10.13182/NT73-A31250},
    }


Dictionary entries (mandatory)
    numberOfPins        37;
    pinDiameter         0.006000;
    pinPitch            0.007600;
    wireDiameter        0.0016;
    wireLeadLen         0.1;        //- I.e. pitch of the wrapper wire
    wetWrapPerimeter    0.175237;   //- Wetted wrapper perimeter (wrapper only)

SourceFiles
    Rehme.C

\*---------------------------------------------------------------------------*/

#ifndef Rehme_H
#define Rehme_H

#include "dragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace dragModels
{

/*---------------------------------------------------------------------------*\
                           Class Rehme Declaration
\*---------------------------------------------------------------------------*/

class Rehme
:
    public dragModel
{
private:

    //- Number of pins in assembly
    label Np_;

    //- Pin diameter
    scalar Dp_;

    //- Pin pitch
    scalar Pp_;

    //- Wrapper wire diameter
    scalar Dw_;

    //- Wrapper wire lead lenght
    scalar Lw_;

    //- Wetted wrapper perimeter
    scalar wetWrapPer_;

    //- Wetted pin perimeter
    scalar wetPinPer_;

    //- Pre-computed coeffs
    scalar A_;
    scalar B1_;
    scalar B2_;

public:

    TypeName("Rehme");

    //- Constructors

        //- Construct for a fluid-structure pair
        Rehme
        (
            const objectRegistry& objReg,
            const dictionary& dict,
            const FSPair& FSPair
        );

        //- Destructor
        virtual ~Rehme(){}


    //- Member Functions

        //- Drag coefficient
        virtual void correctKd(volTensorField& Kd) const;

        //- Scalar drag factor fd in Kd = 0.5*fd*alpha*rho*magU/Dh
        scalar fd(const scalar& Re) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dragModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //