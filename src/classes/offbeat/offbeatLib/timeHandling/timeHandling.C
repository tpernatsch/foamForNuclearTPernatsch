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

\*---------------------------------------------------------------------------*/

// this comment was added to test some git capabilities

#include "timeHandling.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeHandling::timeHandling
(
    const fvMesh& mesh,
    Time& runTime
)
:
    mesh_(mesh),
    runTime_(runTime),
    adjustableTimeStep_(false),
    maxDeltaT_(GREAT),
    minDeltaT_(0),
    lastTimeMarker_(GREAT),
    maxDeltaTMult_(GREAT),
    minDeltaTMult_(0.0),
    adjustableWriteTimeStep_(false),
    writingTimeList_()
{   
    const dictionary& dict(mesh_.time().controlDict());   

    // Time handling variables
    adjustableTimeStep_ = 
        dict.lookupOrDefault<bool>("adjustableTimeStep", false);
    maxDeltaT_ = runTime_.userTimeToTime(
        dict.lookupOrDefault<scalar>("maxDeltaT", GREAT)
    );
    minDeltaT_ = runTime_.userTimeToTime(
        dict.lookupOrDefault<scalar>("minDeltaT", 1)
    );
    maxDeltaTMult_ = 
        dict.lookupOrDefault<scalar>("maxRelativeDeltaTIncrease", GREAT); 
    minDeltaTMult_ = 
        dict.lookupOrDefault<scalar>("minRelativeDeltaTDecrease", 0); 

    // Time writing variables
    adjustableWriteTimeStep_ = 
        dict.lookupOrDefault<bool>("adjustableWriteTimeStep", false);

    if(adjustableWriteTimeStep_)
    {
#ifdef OPENFOAMFOUNDATION
        writingTimeList_.set
        ( 
            new scalarField
            (
                dict.lookup("writeTimeStepList")
            )
        );
#elif OPENFOAMESI
        writingTimeList_.reset
        ( 
            new scalarField
            (
                dict.lookup("writeTimeStepList")
            )
        );
#endif        
    }    
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::timeHandling::setDeltaT(const scalar deltaT)
{
    // Change deltaT only if adjustableTimeStep_ is activated
    if(adjustableTimeStep_)
    {
        // Save old time step size
        scalar oldDeltaT = runTime_.deltaT().value();

        // Set deltaT based on physics requirement (passed to the function
        // via argument) and on maximum rate of change allowed.
        runTime_.setDeltaT
        (
                min
                (
                    deltaT,
                    maxDeltaT_
                )
        );

        // Store temporary value of deltaT
        scalar tmpCurrentDeltaT = runTime_.deltaT().value();

        runTime_.setDeltaT
        (
            max(
                min
                ( 
                    min
                    (
                        tmpCurrentDeltaT,
                        oldDeltaT*maxDeltaTMult_
                    ),
                    maxDeltaT_
                ),
                minDeltaT_
                )
        );

        // Store deltaT after limiting rate of change
        scalar currentDeltaT = runTime_.deltaT().value();

        // Check if next timestep on write list is to be run
        if( adjustableWriteTimeStep_ )
        {
            const scalar predictedNextTime(runTime_.value() + currentDeltaT);

            forAll( writingTimeList_(), elemI )
            {
                // Convert writingTimeList to seconds
                scalar writingTimeI(runTime_.userTimeToTime(
                        writingTimeList_()[elemI]));

                if 
                (
                    runTime_.value() < writingTimeI
                    and 
                    predictedNextTime >= writingTimeI
                )
                {
                    currentDeltaT = writingTimeI - runTime_.value();

                    runTime_.setDeltaT(currentDeltaT);
                    break;
                }
            }
        }

        // Check if next timestep does not go beyond last time marker
        if(runTime_.value() + currentDeltaT > runTime_.endTime().value())
        {
            runTime_.setDeltaT
            (
                runTime_.endTime().value() - runTime_.value()
            );
        }
    }

    // Print deltaT
    Info<< "deltaT = " <<  runTime_.timeToUserTime(runTime_.deltaT().value())
        << " and next timeStep = " 
        << runTime_.timeToUserTime(runTime_.value() + runTime_.deltaT().value()) 
        << endl;
}


bool Foam::timeHandling::write()
{
    if(adjustableWriteTimeStep_)
    { 
        forAll(writingTimeList_(), elemI)
        {
            // Convert to seconds
            scalar timeMarker = 
            runTime_.userTimeToTime(writingTimeList_()[elemI]);

            if
            ( 
                mag(timeMarker - runTime_.value()) == 0
                or
                mag(runTime_.endTime().value() - runTime_.value()) == 0

            )
            {            
                return true;
            }
        }

        return false;
    }
    
    // If adjustable writing is not activated, return true
    return true;

}

// ************************************************************************* //
