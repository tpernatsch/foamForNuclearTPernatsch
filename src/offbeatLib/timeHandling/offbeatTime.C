/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2019 OpenFOAM Foundation
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

#include "offbeatTime.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(offbeatTime, 0);

}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::offbeatTime::timeAdjustment()
{
    if (!deltaTchanged_)
    {
        deltaT_  = userTimeToTime(deltaT_);
    }

    endTime_ = userTimeToTime(endTime_);
    
    if
    (
#ifdef OPENFOAMFOUNDATION    
        writeControl_ == writeControl::runTime
     || writeControl_ == writeControl::adjustableRunTime
#elif OPENFOAMESI
        writeControl_ == writeControls::wcRunTime
     || writeControl_ == writeControls::wcAdjustableRunTime
#endif     
    )    
    {
        writeInterval_ = userTimeToTime(writeInterval_);
    }  

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::offbeatTime::offbeatTime
(
    const word& name,
    const argList& args
)
:
    Time
    (
        name,
        args
    ),
    conversionFactors_(),
    timeUnit_("seconds")
{
    if(controlDict().found("userTime"))
    {
        const dictionary& userTimeDict(controlDict().subDict("userTime"));

        timeUnit_ = userTimeDict.lookupOrDefault<word>("type", "seconds");
    }

    conversionFactors_.insert("seconds", 1.0);
    conversionFactors_.insert("hours", 3600.0);
    conversionFactors_.insert("days", 86400.0);
    
    if (!conversionFactors_.found(timeUnit_))
    {
        FatalErrorInFunction() << "Requested time unit \"" << timeUnit_
            << "\" not recognised. Available units are:"
            << conversionFactors_.sortedToc()
            << abort(FatalError);
    }

    timeAdjustment();

    startTime_  = userTimeToTime(startTime_);
    value()     = userTimeToTime(value());

    deltaTSave_ = deltaT_;
    deltaT0_    = deltaT_;

    Info<< "Selected time unit:" << this->unit() << nl
        << "Selected end time: " << timeToUserTime(endTime().value())
        << this->unit() << nl
        << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::offbeatTime::readDict()
{
    Time::readDict();
    timeAdjustment();
}


bool Foam::offbeatTime::read()
{
    deltaT_  = timeToUserTime(deltaT_);
    endTime_ = timeToUserTime(endTime_);
    if (Time::read())
    {
        timeAdjustment();
        return true;
    }
    else
    {
        return false;
    }
}


Foam::word Foam::offbeatTime::unit() const
{   
    return word(" ") + timeUnit_;
}


Foam::scalar Foam::offbeatTime::userTimeToTime(const scalar t) const
{
    scalar factor = conversionFactors_[timeUnit_];

    return t*factor;
}


Foam::scalar Foam::offbeatTime::timeToUserTime(const scalar t) const
{
    scalar factor = conversionFactors_[timeUnit_];

    return t/factor;
}


Foam::scalar Foam::offbeatTime::timeToUserTimeRatio() const
{
    scalar factor = conversionFactors_[timeUnit_];

    return 1.0/factor;
}


// ************************************************************************* //
