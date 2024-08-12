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

#include "timeDependentAxialProfile.H"
#include "addToRunTimeSelectionTable.H"
#include "IndirectList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(timeDependentAxialProfile, 0);
    
    addToRunTimeSelectionTable
    (
        fastFlux, 
        timeDependentAxialProfile, 
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void Foam::timeDependentAxialProfile::calcAddressing
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


void Foam::timeDependentAxialProfile::calcMinMaxZ
(
    const dictionary& dict
)
{
    // TODO - Since each point can be referenced by multiple cells, this algorithm
    // will likely count each point multiple times. A more efficient approach
    // is needed.
    const fvMesh& referenceMesh = 
    ( 
        mesh_.foundObject<fvMesh>("referenceMesh")
    ) ?
    mesh_.lookupObject<fvMesh>("referenceMesh")
    : 
    mesh_;
    
    const vectorField& points = referenceMesh.points();
    const labelListList& cellPoints = referenceMesh.cellPoints();
    scalarField z = (points & pinDirection_);
    
    zMax_ = -GREAT;
    zMin_ = GREAT;
    
    forAll(addr_, addrI)
    {
        const label cellI = addr_[addrI];
        const labelList& ptIds = cellPoints[cellI];
        
        forAll(ptIds, i)
        {
            const label ptI = ptIds[i];
            zMin_ = min(zMin_, z[ptI]);
            zMax_ = max(zMax_, z[ptI]);
        }
    }
    reduce(zMin_, minOp<scalar>());
    reduce(zMax_, maxOp<scalar>());
    
    if(debug)
    {
        Info<< "Geometry calculated by timeDependentAxialProfile::calcMinMaxZ" << nl
            << tab << "Pin direction: " << pinDirection_ << nl
            << tab << "Z bounds: " << zMin_ << "-" << zMax_ << endl;
    }
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeDependentAxialProfile::timeDependentAxialProfile
(
    const fvMesh& mesh,
    const materials& materials,
    const dictionary& fastFluxOptDict
)
:
    constantFastFlux(mesh, materials, fastFluxOptDict),
    fluxData_(fastFluxOptDict.lookup("fastFlux")),
    timeData_(fastFluxOptDict.lookup("timePoints")),
    method_
    (
        interpolateTableBase::interpolationMethodNames_
        [
            fastFluxOptDict.lookupOrDefault<word>
            ("timeInterpolationMethod", "linear")
        ]
    ),
    fluxTable_(timeData_, fluxData_, method_),
    pinDirection_(),
    axialModel_(),
    addr_()
{
    const globalOptions& globalOpt
    (mesh_.lookupObject<globalOptions>("globalOptions"));
    pinDirection_ = vector(globalOpt.pinDirection());

    calcAddressing(fastFluxOptDict);
    calcMinMaxZ(fastFluxOptDict);

    axialModel_ = axialProfile::New
    (
        mesh, 
        fastFluxOptDict.subDict("axialProfile"), 
        addr_,
        pinDirection_,
        zMin_,
        zMax_
    );
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::timeDependentAxialProfile::~timeDependentAxialProfile()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::timeDependentAxialProfile::correct()
{
    if(mesh_.topoChanging())
    {
        calcAddressing(fastFluxOptDict_);
        calcMinMaxZ(fastFluxOptDict_);
    }

    if(!mesh_.foundObject<volScalarField>("fastFluxPrevIter"))
    {
        flux_.storePrevIter();
    }

    flux_.prevIter().storePrevIter();        
    flux_.storePrevIter();

    const scalarField& fastFluxPrev(flux_.prevIter());
    const scalarField& fastFluxPrevPrev(flux_.prevIter().prevIter());

    // Correct only once per time step or until flux changes
    if
    ( 
        gMax(mag(fastFluxPrevPrev - fastFluxPrev)) > VSMALL 
        or
        mesh_.time().value() > currentTime_
    )
    {           
        currentTime_ = mesh_.time().value();
        // scalar Dt(mesh_.time().deltaT().value());

        // Interpolate fast flux 
        scalar fastFluxAv = fluxTable_
        (mesh_.time().timeToUserTime(currentTime_));

        // Correct the axial profiles
        axialModel_->correct();

        // Build the resulting fast flux and fluence
        scalarField& fastFluxi = flux_.ref();
        const scalarField& f_z = axialModel_->profile();

        forAll(addr_, addrI)
        {
            const label cellI = addr_[addrI];
            
            fastFluxi[cellI] = 
            fastFluxAv*f_z[addrI];
        }

        flux_.correctBoundaryConditions();
        advanceFluence();
    }
    else
    {
        currentTime_ = mesh_.time().value();
    }

}

// ************************************************************************* //
