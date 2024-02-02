/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2306                                                  |
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
    Foam::latticeMap

Description
    Base class for lattices.

Usage
    Exemple for an hexagonal lattice 
    \verbatim
    xyPosLattice
    {
        origin      (0 0);
        nElements   (5 5);
        pitch       0.6;
        type        hexagonal;
        lattice
        (   
            e e F F F      //   y -->
             e F F F F     // x
              F F F F F    // |
               F F F F e   // V
                F F F e e
        );
    }
    \endverbatim

    Exemple for a square lattice 
    \verbatim
    xyPosLattice
    {
        origin      (0 0);
        nElements   (15 15);
        pitch       0.12;
        type        square;
        lattice
        (   
            e e e e e F F F F F e e e e e  //   y -->
            e e e F F F F F F F F F e e e  // x
            e e F F F F F F F F F F F e e  // |
            e F F F F F F F F F F F F F e  // V
            e F F F F F F F F F F F F F e
            F F F F F F F F F F F F F F F
            F F F F F F F F F F F F F F F
            F F F F F F F F F F F F F F F
            F F F F F F F F F F F F F F F
            F F F F F F F F F F F F F F F
            e F F F F F F F F F F F F F e
            e F F F F F F F F F F F F F e
            e e F F F F F F F F F F F e e
            e e e F F F F F F F F F e e e
            e e e e e F F F F F e e e e e
        );
    }
    \endverbatim

\mainauthor 
    Thomas Guilbaud (thomas.guilbaud@epfl.ch)

SourceFiles
    latticeMap.C

\*---------------------------------------------------------------------------*/

#ifndef latticeMap_H
#define latticeMap_H

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class latticeMap Declaration
\*---------------------------------------------------------------------------*/

class latticeMap
{

protected:

    //- Lattice dictionary
    dictionary dict_;

    //- Lattice type, square or hexagonal
    word type_;

    //- X position of the lattice center
    scalar originX_;

    //- Y position of the lattice center
    scalar originY_;
    
    //- FMU pitch distance
    scalar pitch_;

    //- Rotation angle of the lattice, default 0 rad
    scalar rotationAngle_;
    
    //- Number of elements in the lattice on the X axis
    label nx_;

    //- Number of elements in the lattice on the Y axis
    label ny_;

    //- Lattice map in MCNP/Serpent format
    wordList lattice_;

    //- X-position list for all non-empty elements in the lattice map
    scalarList xPos_;

    //- Y-position list for all non-empty elements in the lattice map
    scalarList yPos_;


public:

    //- Constructors
        //- Disallow default bitwise copy construct
        latticeMap(const latticeMap&);

        //- Disallow default bitwise assignment
        void operator=(const latticeMap&);

        //- Constructor from dictionary
        latticeMap(const dictionary& dict);


    //- Destructor
        virtual ~latticeMap(){}


    //- Member functions

        //- Get the X position list
        scalarList getXpositions() const;

        //- Get the Y position list
        scalarList getYpositions() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
