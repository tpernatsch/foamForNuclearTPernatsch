/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2406                                                  |
|    Copyright 2011-2016 OpenFOAM Foundation, 2017-2022 OpenCFD Ltd.         |
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
    Foam::FSDragCoefficientModels::Rehme

Description
    Rheme correlation for pressure drop

    \f[
        P_t = D_p + 1.0444 D_w \\
        P_{w,pin} = N_p \pi (D_p + D_w) \\
        B = \sqrt{\frac{P_t}{D_p}} 
            + \left(
                7.6 \frac{D_p+D_w}{L_w} \left( \frac{P_t}{D_p} \right)^2
            \right)^{2.16} \\
        \\
        A = \frac{P_{w,pin}}{P_{w,pin} + P_{w,wrap}} \\
        B_1 = 64 \sqrt{B} \\
        B_2 = 0.0816 B^{0.9335} \\
        \\
        f_d = A \left( \frac{B_1}{Re} + \frac{B_2}{Re^{0.133}} \right)
    \f]

    where
    \vartable
        N_p        | Number of pins
        D_p        | Pin diameter
        D_w        | Wire diameter
        L_w        | Wire lead length (i.e pitch of the wrapper wire)
        P_{w,pin}  | Wetted pin perimeter
        P_{w,wrap} | Wetted wrapper perimeter
    \endtable

    Example of mandatory dictionary entries:
    \verbatim
    numberOfPins        37;
    pinDiameter         0.006000;
    wireDiameter        0.0016;
    wireLeadLen         0.1;        //- I.e. pitch of the wrapper wire
    wetWrapPerimeter    0.175237;   //- Wetted wrapper perimeter (wrapper only)
    \endverbatim

    Reference publications:
    Klaus Rehme, "PRESSURE DROP CORRELATIONS FOR FUEL ELEMENT SPACERS", Nuclear 
    Technology, 1973, Vol. 17, 1, pages 15-23, DOI:10.13182/NT73-A31250.

    Author: 
    Stefan Radman (sradman@protonmail.com / stefan.radman@epfl.ch)

SourceFiles
    Rehme.C

\*---------------------------------------------------------------------------*/

#ifndef RehmeFSDragCoefficient_H
#define RehmeFSDragCoefficient_H

#include "FSDragCoefficientModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Forward declarations
class FSPair;

namespace FSDragCoefficientModels
{

/*---------------------------------------------------------------------------*\
                           Class Rehme Declaration
\*---------------------------------------------------------------------------*/

class Rehme
:
    public FSDragCoefficientModel
{
protected:

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
            const FSPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~Rehme(){}

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
