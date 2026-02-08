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

#include "constantHeatSource.H"
#include "zeroGradientFvPatchField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(constantHeatSource, 0);

    addToRunTimeSelectionTable
    (
        heatSource,
        constantHeatSource,
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constantHeatSource::constantHeatSource
(
    const fvMesh& mesh,
    const materials& materials,
    const dictionary& heatSourceDict
)
:
    heatSource(mesh, materials, heatSourceDict),
    Q_
    (
        IOobject
        (
            "Q",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Q", dimPower/dimVolume, 0.0),
        zeroGradientFvPatchField<scalar>::typeName
    ),
    currentTime_(mesh_.time().value()),
    maxRelativePowerIncrease_(GREAT),
    maxRelativePowerDecrease_(GREAT),
    maxAbsolutePowerIncrease_(GREAT)
{
    // Adjustable time step functions
    const dictionary& controlDict(mesh_.time().controlDict());

    bool adjustTime = false;
    if (controlDict.found("adjustableTimeStep"))
    {
        adjustTime = controlDict.lookupOrDefault<bool>(
            "adjustableTimeStep", false);
    }
    else
    {
        adjustTime = controlDict.lookupOrDefault<bool>(
            "adjustTimeStep", false);
    }

    if (adjustTime)
    {
        maxRelativePowerIncrease_ =
            (controlDict.lookupOrDefault<scalar>(
                "maxRelativePowerIncrease", GREAT));

        maxRelativePowerDecrease_ =
            (controlDict.lookupOrDefault<scalar>(
                "maxRelativePowerDecrease", GREAT));

        maxAbsolutePowerIncrease_ =
            (controlDict.lookupOrDefault<scalar>(
                "maxAbsolutePowerIncrease", GREAT));
    }

    // Preparae oldTime pointer
    Q_.oldTime();
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::constantHeatSource::~constantHeatSource()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::constantHeatSource::correct()
{
    currentTime_ = mesh_.time().value();
}


Foam::scalar Foam::constantHeatSource::nextDeltaT()
{
    // Calculate deltaT to next time marker
    scalar deltaTtoMarker = this->nextTimeMarker()
            - mesh_.time().value();

    // Starting deltaT for time steppin algorithm cannot be larger than
    // deltaTtoMarker
    scalar deltaT = min(mesh_.time().deltaT().value(), deltaTtoMarker);

    // TODO: decide if relative change must be use or rate
    // Current limiting rate of change
    scalar currentQ(this->averagePowerDensity());
    scalar predictedNewQ(this->predictNextPowerDensity());

    // Calculate relative rate of change
    scalar deltaQ = (predictedNewQ - currentQ);
    scalar relDeltaQ = deltaQ/max(currentQ,SMALL);

    if( currentQ < SMALL)
    {
        // Current heat source is 0. Relative change does not work
        deltaT = this->nextDeltaTAbsoluteChange(deltaT);
    }
    else
    {
        // Current heat source is not 0
        // Check wheter relDeltaQ is negative (power is decreasing)
        if (relDeltaQ < 0)
        {
            // Scale deltaT to allow a maximum relative decrease equal to
            // maxRelativePowerDecrease
            deltaT = deltaT*maxRelativePowerDecrease_ / max(-relDeltaQ, SMALL);

            Info<< "Maximum deltaT calculated by constantHeatSource model: "
                << mesh_.time().timeToUserTime(deltaT)
                << " with a relative power decrease of " << relDeltaQ
                << endl;
        }
        else // Power is increasing or staying constant
        {
            // Scale deltaT to allow a maximum relative increase equal to
            // maxRelativePowerIncrease
            deltaT = deltaT*maxRelativePowerIncrease_ / max(relDeltaQ, SMALL);

            Info<< "Maximum deltaT calculated by constantHeatSource model: "
                << mesh_.time().timeToUserTime(deltaT)
                << " with a relative power increase of " << relDeltaQ
                << endl;
        }
    }

    if (deltaT > deltaTtoMarker)
    {
        deltaT = deltaTtoMarker;

        Info<< "Maximum deltaT re-calculated by constantHeatSource model: "
            << mesh_.time().timeToUserTime(deltaT)
            << " so that deltaT is equal to the deltaT to next time point in heat source history"
            << endl;
    }
    else if (deltaTtoMarker / deltaT < 5)
    {
        // We are approaching the next time marker and adapt the
        // time step to avoid large changes in time step as we approach.
        deltaT = deltaTtoMarker / ceil(deltaTtoMarker/deltaT);

        Info<< "Maximum deltaT re-calculated by constantHeatSource model: "
            << mesh_.time().timeToUserTime(deltaT)
            << " as: deltaTtoMarker < 5*deltaT , i.e. next time point is approaching in heat source history, so time-steps are decreased"
            << endl;
    }

    return deltaT;
}


Foam::scalar Foam::constantHeatSource::nextDeltaTAbsoluteChange(scalar deltaT) const
{
    scalar currentQ(this->averagePowerDensity());
    scalar predictedNewQ(this->predictNextPowerDensity());

    // Calculate absolute rate of change
    scalar deltaQ = (predictedNewQ - currentQ);

    // Set next deltaT
    scalar nextDeltaT = deltaT*maxAbsolutePowerIncrease_/max(deltaQ, SMALL);

    Info<< "Maximum deltaT calculated by constantHeatSource model: "
        << mesh_.time().timeToUserTime(nextDeltaT)
        << " with current heat source " << currentQ <<" W/m3, and an absolute change of " << deltaQ << " W/m3"
        << endl;

    return nextDeltaT;
}

// ************************************************************************* //
