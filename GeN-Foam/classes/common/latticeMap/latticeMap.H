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
    \verbatim
    xyPosLattice
    {
        origin      (0 0);
        nElements   (5 5);
        pitch       0.6;
        type        hexagonal; // square; // hexagonal;
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

    //- Dictionary name
    // word dictName_;

    //- Dictionary
    dictionary dict_;

    //-
    word type_;

    //- 
    scalar originX_;

    //- 
    scalar originY_;
    
    //- 
    scalar pitch_;

    //-
    scalar rotationAngle_;
    
    //- 
    label nx_;

    //-
    label ny_;

    //-
    wordList lattice_;

    //-
    scalarList xPos_;

    //-
    scalarList yPos_;


public:

    //- Constructors
        //- Disallow default bitwise copy construct
        latticeMap(const latticeMap&);

        //- Disallow default bitwise assignment
        void operator=(const latticeMap&);

        //- Constructor from dictionary
        latticeMap(const dictionary& dict);

        //- Constructor from IOobject and subDictname
        latticeMap(IOdictionary object, word latticeMapName);

    //- Destructor
        virtual ~latticeMap(){}


    //- Member functions

        //-
        scalarList getXpositions() const;

        //-
        scalarList getYpositions() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
