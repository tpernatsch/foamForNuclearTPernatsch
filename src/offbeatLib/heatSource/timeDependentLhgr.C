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

#include "timeDependentLhgr.H"
#include "addToRunTimeSelectionTable.H"
#include "IndirectList.H"
    
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(timeDependentLhgr, 0);
    
    addToRunTimeSelectionTable
    (
        heatSource, 
        timeDependentLhgr, 
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeDependentLhgr::timeDependentLhgr
(
    const fvMesh& mesh,
    const materials& materials,
    const dictionary& heatSourceDict
)
:
    constantLhgr(mesh, materials, heatSourceDict, true),
    lhgrData_(heatSourceDict.lookup("lhgr")),
    timeData_(heatSourceDict.lookup("timePoints")),
    method_
    (
        interpolateTableBase::interpolationMethodNames_
        [
            heatSourceDict.lookupOrDefault<word>
            ("timeInterpolationMethod", "linear")
        ]
    ),
    lhgrTable_(timeData_, lhgrData_, method_),
    maxAbsoluteLhgrIncrease_(GREAT)
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
        maxAbsoluteLhgrIncrease_ = 
            (controlDict.lookupOrDefault<scalar>(
                "maxAbsoluteLhgrIncrease", GREAT));
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::timeDependentLhgr::~timeDependentLhgr()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::timeDependentLhgr::correct()
{
    const globalOptions& globalOpt
    (mesh_.lookupObject<globalOptions>("globalOptions"));

    scalar angularFraction = globalOpt.angularFraction();

    if(mesh_.topoChanging())
    {
        calcAddressing(heatSourceDict_);
        calcReferenceDimensions(heatSourceDict_);
    }

    if(!mesh_.foundObject<volScalarField>("QPrevIter"))
    {
        Q_.storePrevIter();
    }

    Q_.prevIter().storePrevIter();        
    Q_.storePrevIter();

    const scalarField& Qprev(Q_.prevIter());
    const scalarField& QprevPrev(Q_.prevIter().prevIter());

    // Correct heat source only once per time step or until power changes
    if
    ( 
        gMax(mag(QprevPrev - Qprev)) > 1e-3 
        or
        mesh_.time().value() > currentTime_
        or
        mesh_.foundObject<volScalarField>("porosity")
    )
    {
        currentTime_ = mesh_.time().value();

        scalarField Vi(mesh_.V(), addr_);
        scalar refVolume(gSum(Vi)/max(refVolumeFraction_, SMALL)); 

        // Interpolate LHGR and calculate average power density Qavg
        scalar lhgr = lookupLHGR(currentTime_);
        scalar Qavg = lhgr*(zMax_ - zMin_)*angularFraction/refVolume;

        // Correct the radial and axial profiles
        radialModel_->correct();
        azimuthalModel_->correct();
        axialModel_->correct();

        // Build the resulting heat source
        scalarField& Qi = Q_.ref();
        scalarField& lhgri = lhgr_.ref();
        const scalarField& f_r = radialModel_->profile();
        const scalarField& h_theta = azimuthalModel_->profile();
        const scalarField& g_z = axialModel_->profile();
        const scalarField& scalingFactors_z = 
            axialModel_->scalingFactorProfile();

        forAll(addr_, addrI)
        {
            const label cellI = addr_[addrI];
            
            Qi[cellI] = Qavg*f_r[addrI]*g_z[addrI]*h_theta[addrI];
            Qi[cellI] *= scalingFactors_z[addrI];  
            lhgri[cellI] = lhgr*g_z[addrI];
        }

        // Reshape power profile due to porosity redistribution
        reshapePowerPorosity();
        
        scalarField Qaddr(Q_.ref(), addr_);        
        scalar powerTot   = lhgr*(zMax_ - zMin_)*angularFraction*gSum(g_z * Vi)/refVolume;
        scalar powerTotOF = gSum( Vi*Qaddr );

        scalar error = 
        100*
        (
            mag(powerTotOF - powerTot)/max(powerTot, VSMALL)
        );

        if(error > 0.01)
        {            
            WarningIn("Foam::timeDependentLhgr::correct()") << nl  
                << "    Total power does not correspond to LHGR*height, "<< nl 
                << "    with an error of "  << error << "%" << nl
                << "    Radial or axial profiles might be not normalized." << nl
                << "    The field Q will be scaled to preserve total power."
                <<  nl << endl;

            radialModel_ -> checkNormalization();
            azimuthalModel_ -> checkNormalization();
            axialModel_ -> checkNormalization();

            scalar scalingFactor(powerTot/powerTotOF);

            forAll(addr_, addrI)
            {
                const label cellI = addr_[addrI];
                
                Qi[cellI] *= scalingFactor;
            }
        }
    }
    else
    {
        currentTime_ = mesh_.time().value();
    }

    Q_.correctBoundaryConditions();
    lhgr_.correctBoundaryConditions();
}


Foam::scalar Foam::timeDependentLhgr::predictNextPowerDensity() const
{
    const globalOptions& globalOpt
    (mesh_.lookupObject<globalOptions>("globalOptions"));

    scalar angularFraction = globalOpt.angularFraction();
    
    // Assume that deltaT does not change
    scalar predictedNextTime = min
    (
        currentTime_ + mesh_.time().deltaT().value(),
        mesh_.time().userTimeToTime
        (
            lhgrTable_.getNextPoint(mesh_.time().timeToUserTime(currentTime_))
        )
    );

    // Volumes of a subset of cells corresponding to addr_ indexing
    scalarField Vi(mesh_.V(), addr_);
    scalar refVolume(gSum(Vi)/max(refVolumeFraction_, SMALL)); 

    // Interpolate LHGR and calculate average power density nextQavg
    scalar lhgr = lookupLHGR(predictedNextTime);
    scalar nextQavg = lhgr*(zMax_ - zMin_)*angularFraction/refVolume;

    // Build the resulting heat source
    const scalarField& f_r = radialModel_->profile();
    const scalarField& h_theta = azimuthalModel_->profile();
    const scalarField& g_z = axialModel_->profile();
    scalarField nextQi = nextQavg * f_r * g_z * h_theta;

    // if(mesh_.foundObject<fvMesh>("referenceMesh"))
    // {
    //     const fvMesh& refMesh =  mesh_.lookupObject<fvMesh>("referenceMesh");

    //     const scalarField& Vref =  refMesh.V();

    //     forAll(addr_, addrI)
    //     {
    //         const label cellI = addr_[addrI];

    //         nextQi[cellI] = nextQi[cellI]
    //                       * Vref[cellI] / Vi[cellI]; 
    //     }
    // }
    
    scalar powerTot   = lhgr*(zMax_ - zMin_)*angularFraction*gSum(g_z * Vi)/refVolume;
    scalar powerTotOF = gSum( Vi*nextQi );
    scalar error = 
    100*
    (
        mag(powerTotOF - powerTot)/max(powerTot, VSMALL)
    );

    if(error > 0.01)
    {           
        scalar scalingFactor(powerTot/powerTotOF);
        nextQi *= scalingFactor;
    }

    return gSum(Vi*nextQi)/gSum(Vi);

}

Foam::scalar Foam::timeDependentLhgr::nextTimeMarker() const
{
    return min
    (
        min
        (
            min
            (
                radialModel_->nextTimeMarker(), 
                azimuthalModel_->nextTimeMarker()
            ),
            mesh_.time().userTimeToTime
            (
                lhgrTable_.getNextPoint(mesh_.time().timeToUserTime(currentTime_))
            )
        ), 
        axialModel_->nextTimeMarker()
    );    
}


Foam::scalar Foam::timeDependentLhgr::lastTimeMarker() const
{
    return min
    (
        min
        (
            min
            (
                radialModel_->lastTimeMarker(), 
                azimuthalModel_->lastTimeMarker()
            ),
            mesh_.time().userTimeToTime(lhgrTable_.getLastPoint())
        ), 
        axialModel_->lastTimeMarker()
    );
}


Foam::scalar Foam::timeDependentLhgr::nextDeltaTAbsoluteChange(scalar deltaT) const
{    
    // Assume that deltaT does not change or that it is limited by lhgr history
    scalar predictedNextTime = min
    (
        currentTime_ + mesh_.time().deltaT().value(),
        mesh_.time().userTimeToTime
        (
            lhgrTable_.getNextPoint(mesh_.time().timeToUserTime(currentTime_))
        )
    );

    scalar currentLhgr(lookupLHGR(currentTime_));
    scalar predictedNewLhgr(lookupLHGR(predictedNextTime));
    
    // Calculate absolute rate of change
    scalar deltaLhgr = (predictedNewLhgr - currentLhgr);

    // Set next deltaT
    scalar nextDeltaT = deltaT*maxAbsoluteLhgrIncrease_/max(deltaLhgr, SMALL);

    Info<< "Maximum deltaT calculated by heatSource model: " 
        << mesh_.time().timeToUserTime(nextDeltaT)
        << " with current lhgr " << currentLhgr <<" W/m, and an absolute change of " << deltaLhgr << " W/m"
        << endl;  

    return nextDeltaT;
}

// ************************************************************************* //
