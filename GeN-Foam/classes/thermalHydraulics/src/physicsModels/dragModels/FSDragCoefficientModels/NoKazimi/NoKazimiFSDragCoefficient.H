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
    Foam::FSDragCoefficientModels::NoKazimi

Description
    NoKazimi correlation for fluid-structure pressure drop
    https://dspace.mit.edu/handle/1721.1/60581

Author
    Stefan Radman (sradman@protonmail.com / stefan.radman@epfl.ch)

Dictionary entries (mandatory)
    pinPitch        0.006;
    pinDiameter     0.0076;
    wireDiameter    0.0016;
    wireLeadLen     0.2;

SourceFiles
    NoKazimi.C

\*---------------------------------------------------------------------------*/

#ifndef NoKazimiFSDragCoefficient_H
#define NoKazimiFSDragCoefficient_H

#include "FSDragCoefficientModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Forward declarations
class FSPair;

namespace FSDragCoefficientModels
{

/*---------------------------------------------------------------------------*\
                           Class NoKazimi Declaration
\*---------------------------------------------------------------------------*/

class NoKazimi
:
    public FSDragCoefficientModel
{
protected:

    //- Pre-computed coeffs
    scalar A_;
    scalar B_;
    scalar C_;

public:

    TypeName("NoKazimi");

    //- Constructors

        //- Construct for a fluid-structure pair
        NoKazimi
        (
            const FSPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~NoKazimi(){}

    //- Member Functions

        //- Drag coefficient
        virtual scalar value(const label& celli) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FSDragCoefficientModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //