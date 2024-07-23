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
    Foam::timeProfile

Description
    Base class for time profile functions. This object embeds the OpenFOAM 
    Function1 object used to represents tables. timeProfile is also able to
    input data from an FMI.

    Example of time dependant source entries:
    \verbatim
        // Using table
        externalReactivityTimeProfile
        {
            type        table;
            startTime   100;
            table       table
            (
                (   0   0)
                (   1   0.0005767131    ) //- i.e. 57.67131 pcm, == 0.2*Beta
            );
        }

        // Using FMI
        externalReactivityTimeProfile
        {
            type            fmi;
            nameFromFMU     gfExtReact;
            initialValue    0.0;
        }
    \endverbatim

    Author of this file (and associated .C or included .H files): 
    Thomas Guilbaud <thomas.guilbaud@epfl.ch;>, EPFL/Transmutex SA (Swizerland)

SourceFiles
    timeProfile.C

\*---------------------------------------------------------------------------*/

#ifndef timeProfile_H
#define timeProfile_H

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class timeProfile Declaration
\*---------------------------------------------------------------------------*/

class timeProfile
{

protected:

    //- Dictionary
    dictionary dict_;

    //- Runtime address
    const Time& runTime_;

    //- Type of the function used for calculation (e.g. table, ...)
    word type_;

    //- Time offset
    scalar startTime_;

    //-
    scalar oldTime_;

    //-
    scalar currentValue_;

    //-
    scalar oldValue_;

    //- Function pointer
    autoPtr<Function1<scalar>> functionPtr_;

    //- Flag true if FMI port correctly setup
    bool isfmiPortSet_;

    //- Name of the FMI port
    word nameFromFMU_;

public:

    //- Constructors
        //- Disallow default bitwise copy construct
        // timeProfile(const timeProfile&);

        //- Disallow default bitwise assignment
        // void operator=(const timeProfile&);

        //- Constructor from dictionary
        timeProfile(const dictionary& dict, const Time& runTime);

        //- Constructor from IOobject and subDictname
        timeProfile(
            IOdictionary object,
            word timeProfileName,
            const Time& runTime
        );

    //- Destructor
        virtual ~timeProfile(){}


    //- Member functions

        //- Initialiye the FMI port
        void initializeFMI();

        //- Test if Function1 has been declared or if FMI port has been create
        bool valid() const;

        //- Get the value at a given time minus the start time offset.
        scalar value(scalar time);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
