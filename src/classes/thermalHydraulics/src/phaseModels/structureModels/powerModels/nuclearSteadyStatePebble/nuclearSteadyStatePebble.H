/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2406                                                  |
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
    Foam::powerModels::nuclearSteadyStatePebble

Description
    Describes triso and pebble of a pebble-bed reactor. Only steady-state.

    Authors:  
    Yvver Robert, Ludovic Jantzen (UC Berkeley)
    Review: Carlo Fiorina
    
    

Usage
    \verbatim
            type                nuclearSteadyStatePebble;
            
            // It can be shown that by averaging a/multiple pebble(s) over
            // a volume of any shape, the interfacialArea and volumeFraction
            // of the resulting porous pebble structure are not independent, yet
            // are tied by volumetricArea = 3*volumeFraction/outerPinRadius
            powerDensityNeutronics            19.96567386e6; 
            pebbleSurfaceTemperatureGuess 939;

            pebbleCoreRadius        1.38e-2;
            pebbleMatrixRadius      1.8e-2;
            pebbleShellRadius       2e-2;
            trisoFuelRadius         212.5e-6;
            trisoBufferRadius       312.5e-6;
            trisoInnerPyCRadius     352.5e-6;
            trisoSiCRadius          387.5e-6;
            trisoOuterPyCRadius     427.5e-6;

            nTRISO                  9022;

            trisoFuelKCoeffs      (3.3073 0 0 0 0 0 0 0);
            trisoBufferKCoeffs    (0.50 0 0 0 0 0 0 0); 
            trisoPyCKCoeffs       (4.00 0 0 0 0 0 0 0);
            trisoSiCKCoeffs       (90.3 0 0 0 0 0 0 0); 
            pebbleGraphiteKCoeffs (41.68 0 0 0 0 0 0 0);

            trisoFuelDensityCoeffs      (11031.432 -0.17907428 -9.1359257e-05 2.2133072e-09 0 0 0 0); thesis
            trisoBufferDensityCoeffs    (1000 0 0 0 0 0 0 0);
            trisoPyCDensityCoeffs       (1900 0 0 0 0 0 0 0);
            trisoSiCDensityCoeffs       (4.015E+02 -5.995E-01 3.898E-04 -9.380E-08 0 0 0 0);
            pebbleGraphiteDensityCoeffs (8.12E+01 -6.14E-02 2.73E-05 -5.18E-09 0 0 0 0);
            coolantDensityCoeffs        (2413.03 -0.4884 0 0 0 0 0 0);

    \endverbatim

SourceFiles
    nuclearSteadyStatePebble.C

\*---------------------------------------------------------------------------*/

#ifndef nuclearSteadyStatePebble_H
#define nuclearSteadyStatePebble_H

#include "powerModel.H"
#include "IOFieldFields.H"
#include "Function1.H"
#include "Polynomial.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace powerModels
{

/*---------------------------------------------------------------------------*\
                           Class nuclearSteadyStatePebble Declaration
\*---------------------------------------------------------------------------*/

class nuclearSteadyStatePebble
:
    public powerModel
{
protected:

    //- Fraction of power from neutronics that goes to this powerModel
    //  This take sinto account for instance cases where some of the 
    //  heating goes to the moderator
    scalarList fractionOfPowerFromNeutronics_;

    //- First guess for pebble surface temperature, to avoid initializing temperatures to 0K
    scalarList TpS_guess_;


    //- Average fuel and matrix temperatures, and pebble, matrix and fuel surface temperatures
    volScalarField TpS_;
    volScalarField Tmav_;
    volScalarField Tmout_;
    volScalarField TfS_;
    volScalarField Tfav_;
    volScalarField Tfmax_;

    //- Scalars that are input in the IOdictionary for passing min/max temperatures
    scalar Tfavmax_;
    scalar Tfavmin_;
    scalar Tfmaxmax_;
    scalar Tfmaxmin_;
    scalar Tmoutmax_;
    scalar Tmoutmin_;
    scalar Tmavmax_;
    scalar Tmavmin_;
    scalar TpSmax_;
    scalar TpSmin_;
    scalar TfSmax_;
    scalar TfSmin_;

    //- Pebble radii
    scalarList rcore_;   // Pebble inner graphite core
    scalarList rmatrix_; // Pebble fuel layer with matrix
    scalarList rshell_;  // Pebble outer graphite shell

    //- TRISO radii
    scalarList rfuel_;   // TRISO fuel kernel
    scalarList rbuffer_; // TRISO buffer layer
    scalarList rPyCin_;  // TRISO inner PyC layer
    scalarList rSiC_;    // TRISO SiC layer
    scalarList rPyCout_; // TRISO outer PyC layer

    //- Thermal conductivities
    List< Polynomial<8> > kfuelCoeffs_;
    List< Polynomial<8> > kbufferCoeffs_;
    List< Polynomial<8> > kSiCCoeffs_;
    List< Polynomial<8> > kPyCCoeffs_;
    List< Polynomial<8> > kgraphiteCoeffs_;
    scalarList kfuel_;     // TRISO fuel kernel
    scalarList kbuffer_;   // TRISO buffer layer
    scalarList kPyC_;      // TRISO PyC layers
    scalarList kSiC_;      // TRISO SiC layer
    scalarList kgraphite_; // Pebble matrix graphite

    //- Matrix effective thermal conductivity
    volScalarField keffMatrix_;

    //- Densities
    List< Polynomial<8> > rhofuelCoeffs_;
    List< Polynomial<8> > rhobufferCoeffs_;
    List< Polynomial<8> > rhoSiCCoeffs_;
    List< Polynomial<8> > rhoPyCCoeffs_;
    List< Polynomial<8> > rhographiteCoeffs_;
    List< Polynomial<8> > rhocoolantCoeffs_;
    volScalarField rhofuel_;     // TRISO fuel kernel
    volScalarField rhobuffer_;   // TRISO buffer layer
    volScalarField rhoPyC_;      // TRISO PyC layers
    volScalarField rhoSiC_;      // TRISO SiC layer
    volScalarField rhographite_; // Pebble matrix graphite
    volScalarField rhocoolant_; // Coolant

    //- Pebble/fluid surface heat flux
    volScalarField surfaceHeatFlux_;

    //- Number of TRISO per pebble
    labelList Ntriso_;

    //- List that maps from cell index to region index
    labelList cellToRegion_;

    //- List that maps the region index to the region name
    wordList regionIndexToRegionName_;

    //-
    wordList regionNames_;

    //- Just to avoid typing constant::mathematical::pi every time...
    scalar pi_;

public:

    //- Runtime type information
    TypeName("nuclearSteadyStatePebble");

    // Constructors

        nuclearSteadyStatePebble
        (
            structure& structure,
            const dictionary& dicts
        );

        //- Destructor
        virtual ~nuclearSteadyStatePebble();

    // Member Functions

        //- Set interfacial area consistently with the pin geometry. In
        //  particular, if the overall total volume-fraction of the pins
        //  in a given mesh cell is alpha, it can be shown that the interfacial
        //  area density iA of the pin surface is 2*alpha/outerRadius
        virtual void setInterfacialArea();

        void updateLocalAvgGlobalMinMaxT
        (
            scalar& T,
            scalar& Tmin,
            scalar& Tmax
        );

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
