/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2020 OpenFOAM Foundation
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
    Foam::offbeatTime

Description
    Manage time in OFFBEAT in days/hours/seconds as decided by the user in 
    controlDict.

    The user should be aware that the choice of time unit must be consistent 
    across the input files, i.e. all time dependent tables provided in the 
    solverDict or in the boundary conditions must be provided accordingly to the 
    time unit selected in the controlDict. 

Usage
    In controlDict    
    \verbatim    
    ...

    userTime
    {
        // type is seconds by default
        type seconds;//hours;//days;
    }

    ...
    \endverbatim  

SourceFiles
    offbeatTime.C

\*---------------------------------------------------------------------------*/

#ifndef offbeatTime_H
#define offbeatTime_H

#include "Time.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class offbeatTime Declaration
\*---------------------------------------------------------------------------*/

class offbeatTime
:
    public Time
{
    // Private data

    // Private Member Functions

        //- Adjust read time values
        void timeAdjustment();

public:

    //- Runtime type information
    TypeName("offbeatTime");

    // Public data
    
        //- Conversion factors from user time to time
        HashTable<scalar, word> conversionFactors_;

        //- User-selected time units
        word timeUnit_;


    // Constructors

        //- Construct from objectRegistry arguments
        offbeatTime
        (
            const word& name,
            const argList& args
        );

        //- Disallow default bitwise copy construction
        offbeatTime(const offbeatTime&) = delete;


    //- Destructor
    virtual ~offbeatTime()
    {}


    // Member Functions

        // Access

            //- Return time unit
            virtual word unit() const;

            //- Current user time
            inline scalar userTime() const
            {
                return timeToUserTime(this->value());
            }
            
        // Member Functions overriding the virtual functions in time

            //- Convert the user-time (e.g. hours) to real-time (s).
            virtual scalar userTimeToTime(const scalar theta) const;

            //- Convert the real-time (s) into user-time (e.g. hours)
            virtual scalar timeToUserTime(const scalar t) const;

            //- Ratio between real-time and user-time
            virtual scalar timeToUserTimeRatio() const;

            //- Read the control dictionary and set the write controls etc.
            virtual void readDict();


        // Edit

            //- Read the controlDict and set all the parameters
            virtual bool read();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const offbeatTime&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
