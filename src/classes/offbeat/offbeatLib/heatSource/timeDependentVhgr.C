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

#include "timeDependentVhgr.H"
#include "addToRunTimeSelectionTable.H"
#include "IndirectList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(timeDependentVhgr, 0);

    addToRunTimeSelectionTable
    (
        heatSource,
        timeDependentVhgr,
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void Foam::timeDependentVhgr::calcAddressing
(
    const dictionary& dict
)
{
    addr_ = labelList();

    wordList zoneNames(dict.lookup("materials"));

    forAll(zoneNames, i)
    {
        const label zoneI = mesh_.cellZones().findZoneID(zoneNames[i]);

        if (zoneI == -1)
        {
            FatalIOErrorInFunction(dict.lookup("materials"))
                << "cellZone " << zoneNames[i] << " not found on mesh."
                << abort(FatalIOError);
        }

        addr_.append(mesh_.cellZones()[zoneI]);
    }
}



// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeDependentVhgr::timeDependentVhgr
(
    const fvMesh& mesh,
    const materials& materials,
    const dictionary& heatSourceOptDict
)
:
    heatSource(mesh, materials, heatSourceOptDict),
    vhgrData_(heatSourceOptDict.lookup("vhgr")),
    timeData_(heatSourceOptDict.lookup("timePoints")),
    method_
    (
        interpolateTableBase::interpolationMethodNames_
        [
            heatSourceOptDict.lookupOrDefault<word>
            ("timeInterpolationMethod", "linear")
        ]
    ),
    vhgrTable_(timeData_, vhgrData_, method_),
    addr_()
{
    calcAddressing(heatSourceOptDict);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::timeDependentVhgr::~timeDependentVhgr()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::timeDependentVhgr::correct()
{
    if(mesh_.topoChanging())
    {
        calcAddressing(heatSourceOptDict_);
    }

    if(!mesh_.foundObject<volScalarField>("QPrevIter"))
    {
        Q_.storePrevIter();
    }

    Q_.prevIter().storePrevIter();
    Q_.storePrevIter();

    const scalarField& Qprev(Q_.prevIter());
    const scalarField& QprevPrev(Q_.prevIter().prevIter());

    //-Correct heat source only once per time step or until power changes
    if
    (
        gMax(mag(QprevPrev - Qprev)) > 1e-3
        or
        mesh_.time().value() > currentTime_
    )
    {
        currentTime_ = mesh_.time().value();

        scalarField Vi(mesh_.V(), addr_);

        // Interpolate vhgr and calculate average power density Qavg
        scalar vhgr = lookupVHGR(currentTime_);
        scalar Qavg = vhgr;

        // Build the resulting heat source
        scalarField& Qi = Q_.ref();
        forAll(addr_, addrI)
        {
            const label cellI = addr_[addrI];
            Qi[cellI] = Qavg;
        }
    }
    else
    {
        currentTime_ = mesh_.time().value();
    }

    Q_.correctBoundaryConditions();
}


Foam::scalar Foam::timeDependentVhgr::predictNextPowerDensity() const
{
    //- Assume that deltaT does not change
    scalar predictedNextTime = min
    (
        currentTime_ + mesh_.time().deltaT().value(),
        mesh_.time().userTimeToTime
        (
            vhgrTable_.getNextPoint(mesh_.time().timeToUserTime(currentTime_))
        )
    );

    // Volumes of a subset of cells corresponding to addr_ indexing
    scalarField Vi(mesh_.V(), addr_);

    // Interpolate vhgr and calculate average power density nextQavg
    scalar vhgr = lookupVHGR(predictedNextTime);
    scalar nextQavg = vhgr;

    return nextQavg;

}


Foam::scalar Foam::timeDependentVhgr::nextDeltaTAbsoluteChange(scalar deltaT) const
{    
    // Assume that deltaT does not change or that it is limited by vhgr history
    scalar predictedNextTime = min
    (
        currentTime_ + mesh_.time().deltaT().value(),
        mesh_.time().userTimeToTime
        (
            vhgrTable_.getNextPoint(mesh_.time().timeToUserTime(currentTime_))
        )
    );

    scalar currentVhgr(lookupVHGR(currentTime_));
    scalar predictedNewVhgr(lookupVHGR(predictedNextTime));
    
    // Calculate absolute rate of change
    scalar deltaVhgr = (predictedNewVhgr - currentVhgr);

    // Set next deltaT
    scalar nextDeltaT = deltaT*maxAbsolutePowerIncrease_/max(deltaVhgr, SMALL);

    Info<< "Maximum deltaT calculated by heatSource model: " 
        << mesh_.time().timeToUserTime(nextDeltaT)
        << " with current vhgr " << currentVhgr <<" W/m3, and an absolute change of " << deltaVhgr << " W/m3"
        << endl;  

    return nextDeltaT;
}

Foam::scalar Foam::timeDependentVhgr::averagePowerDensity() const
{
    scalarField Vi(mesh_.V(), addr_);
    scalarField Qi(Q_, addr_);

    return gSum(Vi*Qi)/gSum(Vi);
}

Foam::scalar Foam::timeDependentVhgr::nextTimeMarker() const
{
    return vhgrTable_.getNextPoint(mesh_.time().timeToUserTime(currentTime_));
}


Foam::scalar Foam::timeDependentVhgr::lastTimeMarker() const
{
    return mesh_.time().userTimeToTime(vhgrTable_.getLastPoint());
}

// ************************************************************************* //
