/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2412                                                  |
|    Copyright 2011-2016 OpenFOAM Foundation, 2017-2024 OpenCFD Ltd.          |
-------------------------------------------------------------------------------
License
    This file is part of GeN-Foam.

    GeN-Foam is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    GeN-Foam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    This offering is not approved or endorsed by the OpenFOAM Foundation nor
    OpenCFD Limited, producer and distributor of the OpenFOAM(R)software via
    www.openfoam.com, and owner of the OPENFOAM(R) and OpenCFD(R) trademarks.

    This particular snippet of code is developed according to the developer's
    knowledge and experience in OpenFOAM. The users should be aware that
    there is a chance of bugs in the code, though we've thoroughly test it.
    The source code may not be in the OpenFOAM coding style, and it might not
    be making use of inheritance of classes to full extent.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::powerModels::interpolatedNuclearFuelPin

Description
    Model for representing a nuclear fuel pin with constant material
    properties. Models fuel, gap, cladding and solves for the heat equation
    with a finite differences scheme. Only radial heat diffusion is considered.
    Derivation of the equations is provided at the end of the present file.

    Authors:
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com),
    Carlo Fiorina (carlo.fiorina@outlook.com / carlo.fiorina@epfl.ch)

SourceFiles
    interpolatedNuclearFuelPin.C

\*---------------------------------------------------------------------------*/

#ifndef interpolatedNuclearFuelPin_H
#define interpolatedNuclearFuelPin_H

#include "powerModel.H"
#include "IOFieldFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace powerModels
{

/*---------------------------------------------------------------------------*\
                           Class interpolatedNuclearFuelPin Declaration
\*---------------------------------------------------------------------------*/

class interpolatedNuclearFuelPin
:
    public powerModel
{
protected:

    //- Fuel outer radius
    scalarList rfo_;

    //- List that maps from cell index to region index
    labelList cellToRegion_;

    //- List that maps the region index to the region name
    wordList regionIndexToRegionName_;

    //-
    wordList regionNames_;

    volScalarField interpolatedHeatFlux_;
    volScalarField interpolatedPowerDensity_;

    //- Structure temperature
    volScalarField T_;


public:

    //- Runtime type information
    TypeName("interpolatedNuclearFuelPin");

    // Constructors

        interpolatedNuclearFuelPin
        (
            structure& structure,
            const dictionary& dicts
        );

        //- Destructor
        virtual ~interpolatedNuclearFuelPin();

    // Member Functions

        //- Set interfacial area consistently with the pin geometry. In
        //  particular, if the overall total volume-fraction of the pins
        //  in a given mesh cell is alpha, it can be shown that the interfacial
        //  area density iA of the pin surface is 2*alpha/outerRadius
        virtual void setInterfacialArea();

        //- In two-phase solvers :
        //  HT = H1*frac1*Tfluid1 + H2*frac2*Tfluid2
        //  H =  H1*frac1 + H2*frac2
        //- In a mono-phase solver :
        //  HT = H*Tfluid
        //  H = H
        virtual void correct
        (
            const volScalarField& HT,
            const volScalarField& H
        );

        //- Return surface temperature (i.e. outer cladding temperature)
        virtual void correctT(volScalarField& T) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace thermoSubstructures
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //