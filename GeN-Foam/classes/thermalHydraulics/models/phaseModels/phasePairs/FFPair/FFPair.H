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

Description
    Class to handle dimensionless numbers that characterise the system, e.g.
    Reynolds, etc. While some can (rightfully) think that some of these numbers
    should be directly members of other classes (e.g., Reynolds of fluids
    should be in the fluid class), in general these numbers depend on a
    multitude of parameters. You cannot define a Reynolds for a fluid by
    itself, it needs to be geometrically constratined in a certain domain 
    (so that a hydraulic diameter can be defined), so it is not properly
    "just" of the fluid. While this logic can of course be debate, I wanted
    to keep it simple and not start dispersing dimensionless numbers randomly
    around the code. 

Author
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

Class
    Foam::FFPair

SourceFiles
    FFPair.C

\*---------------------------------------------------------------------------*/

#ifndef FFPair_H
#define FFPair_H

#include "phasePair.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Forward declarations
class fluid;

/*---------------------------------------------------------------------------*\
                           Class FFPair Declaration
\*---------------------------------------------------------------------------*/

class FFPair
    : public phasePair
{
protected :

    //- Ref to fluid1
    const fluid& fluid1_;

    //- Ref to fluid2
    const fluid& fluid2_;

    //- Phase fraction of the dispersed phase
    volScalarField alphaDispersed_;

    //- Hydraulic diameter of the dispersed phase
    volScalarField DhDispersed_;

    //- Hydraulic diameter of the continuous phase
    volScalarField DhContinuous_;

    //- Density of the continuous phase
    volScalarField rhoContinuous_;

    //- Effective momentum diffusivity of the continuous phase
    volScalarField nuContinuous_;

    //- Pr of the continuous phase
    volScalarField PrContinuous_;

    //- Reynolds of the pair
    volScalarField Re_;

    //- Minimum allowable Reynolds
    dimensionedScalar minRe_;

    //- Magnitude of relative velocity difference between phases
    volScalarField magUr_;

public:

    TypeName("FFPair");

    //- Constructors

        FFPair
        (
            const fluid& fluid1,
            const fluid& fluid2
        );

        //- Destructor
        virtual ~FFPair(){}

    //- Member Functions

        //- Update dimensionless numbers
        void correct();

        //- Access

            const fluid& fluid1() const
            {
                return fluid1_;
            }

            const fluid& fluid2() const
            {
                return fluid2_;
            }

            const volScalarField& alphaDispersed() const
            {
                return alphaDispersed_;
            }

            const volScalarField& DhDispersed() const
            {
                return DhDispersed_;
            }

            const volScalarField& DhContinuous() const
            {
                return DhContinuous_;
            }

            const volScalarField& rhoContinuous() const
            {
                return rhoContinuous_;
            }

            const volScalarField& nuContinuous() const
            {
                return nuContinuous_;
            }

            const volScalarField& PrContinuous() const
            {
                return PrContinuous_;
            }

            const volScalarField& Re() const
            {
                return Re_;
            }       

            const volScalarField& magUr() const
            {
                return magUr_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
