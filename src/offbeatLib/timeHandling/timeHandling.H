/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::timeHandling

Description
    Class that handles variation in time step during simulation.

SourceFiles
    timeHandling.C

Usage
    In controlDict:
    \verbatim

    ...

    adjustableTimeStep true;

    maxDeltaT         1e6;
    minDeltaT         1e-6;

    maxRelativeDeltaTIncrease 1e9;
    minRelativeDeltaTDecrease 1e9;

    maxRelativePowerIncrease  0.25;
    maxRelativePowerDecrease  0.25;

    maxBurnupIncrease 0.1;
    maxAverageCreep   1e-4;
    maxMaximumCreep   1e-4;
    maxFGR            1e9;

    // Allow adjustable time-step writing or not 
    adjustableWriteTimeStep         off;
    // writeTimeStepList
    // (
    //  0.0
    //  91822464 
    //  92822464 
    //  100462464
    // );

    \endverbatim

\todo
    Move all adjustable time options in solverDict

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)\n

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef timeHandling_H
#define timeHandling_H

#include "offbeatTime.H"
#include "fvMesh.H"
#include "IOdictionary.H"
#include "Table.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class timeHandling Declaration
\*---------------------------------------------------------------------------*/

class timeHandling
{
    // Private data
    
        //- Reference to object fvMesh
        const fvMesh& mesh_;

        //- Non-const reference to runTime
        offbeatTime& runTime_;

        //- Adjustable time-step 
        bool adjustTimeStep_;

        //- Maximum time increment
        scalar maxDeltaT_;

        //- Minimum time increment
        scalar minDeltaT_;

        //- Last time marker;
        scalar lastTimeMarker_;

        //- Multiplier for maximum time step increase
        scalar maxDeltaTMult_;

        // //- Multiplier for maximum time step decrease
        // scalar minDeltaTMult_;

        //- Adjustable time-step writing
        bool adjustWriteTimeStep_;

        //- List of writing time steps 
        autoPtr<scalarField> writingTimeList_;
        
    // Private Member Functions

        //- Disallow default bitwise copy construct
        timeHandling(const timeHandling&);

        //- Disallow default bitwise assignment
        void operator=(const timeHandling&);

protected:
    
    // Protected data  
        
public:

    // Constructors

        //- Construct from mesh
        timeHandling(const fvMesh&, offbeatTime& runTime);
        

    //- Destructor
    ~timeHandling()
    {};


    // Access functions


    // Member Functions

        //- Set time step size (in seconds)
        void setDeltaT(const scalar deltaT);

        //- Set last time marker (in seconds) after which the simulation is stopped
        void setLastTimeMarker(scalar timeMarker)
        {
            lastTimeMarker_ = timeMarker;
        }

        //- Return true if a write is pending
        bool write();
    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
