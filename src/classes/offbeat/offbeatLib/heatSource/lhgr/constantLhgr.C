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

#include "constantLhgr.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchField.H"
#include "IndirectList.H"
#include "fvc.H"
    
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(constantLhgr, 0);
    
    addToRunTimeSelectionTable
    (
        heatSource, 
        constantLhgr, 
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::constantLhgr::calcAddressing
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


void Foam::constantLhgr::calcReferenceDimensions
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
    scalarField r = mag(points - (z*pinDirection_));

    zMax_ = -GREAT;
    zMin_ = GREAT;
    rMax_ = 0;
    
    forAll(addr_, addrI)
    {
        const label cellI = addr_[addrI];
        const labelList& ptIds = cellPoints[cellI];
        
        forAll(ptIds, i)
        {
            const label ptI = ptIds[i];
            zMin_ = min(zMin_, z[ptI]);
            zMax_ = max(zMax_, z[ptI]);
            rMax_ = max(rMax_, r[ptI]);
        }
    }
    reduce(zMin_, minOp<scalar>());
    reduce(zMax_, maxOp<scalar>());
    reduce(rMax_, maxOp<scalar>());
    
    if(debug)
    {
        Info<< "Geometry calculated by timeDependentLhgr::calcReferenceDimensions" << nl
            << tab << "Pin direction: " << pinDirection_ << nl
            << tab << "Z bounds: " << zMin_ << "-" << zMax_ << nl
            << tab << "Radius: "   << rMax_ << nl << endl;
    }
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constantLhgr::constantLhgr
(
    const fvMesh& mesh,
    const materials& materials,
    const dictionary& heatSourceOptDict,
    const bool timeDependent
)
:
    heatSource(mesh, materials, heatSourceOptDict),
    lhgr_
    (
        IOobject
        (
            "lhgr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("lhgr", dimPower/dimLength, 0.0),
        zeroGradientFvPatchField<scalar>::typeName
    ),
    constLhgr_(),
    pinDirection_(),
    radialDirection_(),
    radialModel_(),
    azimuthalModel_(),
    axialModel_(),
    addr_(),
    updateFromPK_(heatSourceOptDict.getOrDefault<bool>("updateFromPK", false)),
    powerOld_(heatSourceOptDict.getOrDefault<scalar>("initialPower", 0))

{
    if( !timeDependent )
    {
        heatSourceOptDict.lookup("lhgr") >> constLhgr_;
    }

    const globalOptions& globalOpt
    (mesh_.lookupObject<globalOptions>("globalOptions"));
    pinDirection_ = vector(globalOpt.pinDirection());
    radialDirection_ = vector(globalOpt.radialDirection());

    calcAddressing(heatSourceOptDict);
    calcReferenceDimensions(heatSourceOptDict);
    
    radialModel_ = radialProfile::New
    (
        mesh, 
        heatSourceOptDict.subDict("radialProfile"), 
        addr_,
        pinDirection_,
        rMax_
    );
    
    azimuthalModel_ = azimuthalProfile::New
    (
        mesh, 
        heatSourceOptDict.subOrEmptyDict("azimuthalProfile"), 
        addr_,
        pinDirection_,
        radialDirection_
    );

    axialModel_ = axialProfile::New
    (
        mesh, 
        heatSourceOptDict.subDict("axialProfile"), 
        addr_,
        pinDirection_,
        zMin_,
        zMax_
    );
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::constantLhgr::~constantLhgr()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::constantLhgr::correct()
{
    const globalOptions& globalOpt
    (mesh_.lookupObject<globalOptions>("globalOptions"));

    scalar angularFraction = globalOpt.angularFraction();

    if(mesh_.topoChanging())
    {
        calcAddressing(heatSourceOptDict_);
        calcReferenceDimensions(heatSourceOptDict_);
    }

    if(!mesh_.foundObject<volScalarField>("QPrevIter"))
    {
        Q_.storePrevIter();
    }

    Q_.prevIter().storePrevIter();        
    Q_.storePrevIter();

    // Correct heat source only once per time step or until power changes
    if
    ( 
        mesh_.time().value() > currentTime_
        or
        mesh_.foundObject<volScalarField>("porosity")
        or updateFromPK_
    )
    {  
        currentTime_ = mesh_.time().value();

        scalarField Vi(mesh_.V(), addr_);
        scalar totalVol(gSum(Vi)/angularFraction);

        // Interpolate LHGR and calculate average power density Qavg
        scalar Qavg = constLhgr_*(zMax_ - zMin_)/totalVol;
        
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


        //- If the simulation is run coupled with pointKinetics, it could be useful to update the powerShape to account for the reactivity effects
        //- This is done already in GeN-Foam. If two power fields exist, i.e. powerDensityOld and powerDensityNew, then if i compute 
        // - fvc::integrate(powerDensityNew)/fvc::integrate(powerDenistyOld) this gives me the PK scaling factor

        if (updateFromPK_)
        {
            const word PKRegion(heatSourceOptDict_.get<word>("PKRegionName"));
            const fvMesh& PKMesh(mesh_.time().lookupObject<fvMesh>(PKRegion));
            const scalar powerNew(fvc::domainIntegrate(PKMesh.lookupObject<volScalarField>("powerDensityNeutronics")).value());
            if(powerOld_!=0)
            {
                Q_*=powerNew/powerOld_;
                lhgr_*=powerNew/powerOld_;
                constLhgr_*=powerNew/powerOld_;
                Info << nl<<"Power changed from "<< powerOld_ <<" W to " << powerNew << " W"<<endl;
                Info <<"Rescaling the power by a factor " << powerNew/powerOld_<<nl<<endl;
            }
            powerOld_=powerNew;
        }

        forAll(addr_, addrI)
        {
            const label cellI = addr_[addrI];
            
            Qi[cellI] = Qavg*f_r[addrI]*g_z[addrI]*h_theta[addrI];
            lhgri[cellI] = constLhgr_*g_z[addrI];
        }

        // Reshape power profile due to porosity redistribution
        reshapePowerPorosity();

        scalarField Qaddr(Q_.ref(), addr_);
        scalar powerTot   = constLhgr_*(zMax_ - zMin_)*angularFraction;
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


Foam::scalar Foam::constantLhgr::predictNextPowerDensity() const
{
    return averagePowerDensity();
}


Foam::scalar Foam::constantLhgr::averagePowerDensity() const
{
    scalarField Vi(mesh_.V(), addr_);
    scalarField Qi(Q_, addr_);

    return gSum(Vi*Qi)/gSum(Vi);
}

void Foam::constantLhgr::reshapePowerPorosity()
{
    // Scale volumetric power according to porosity profile
    if
    (
        mesh_.foundObject<volScalarField>("porosity")
    )
    {
        // Ref to Q field
        scalarField& Qi = Q_.ref();
        
        // Initialize avg porosity in the domain
        scalar avgPorosity(0); 

        // Local porosity field
        scalarField porosity(
            mesh_.lookupObject<volScalarField>("porosity").internalField(), addr_);

        // Local volume
        scalarField volume(mesh_.V().field(), addr_);

        // Compute average porosity
        avgPorosity = gSum(porosity*volume)/gSum(volume);

        forAll(addr_, addrI)
        {
            const label cellI = addr_[addrI];

            // Get local value of porosity (bounded in [0,1])
            const scalar pI = max(min(porosity[addrI], 1.0), 0.0);

            // Scale factor for volumetric power
            Qi[cellI] *= (1-pI)/(1-avgPorosity);

        }
    }
}


Foam::scalar Foam::constantLhgr::nextTimeMarker() const
{
    return min
    (
        min
        (
            radialModel_->nextTimeMarker(), 
            azimuthalModel_->nextTimeMarker()
        ), 
        axialModel_->nextTimeMarker()
    );
}


Foam::scalar Foam::constantLhgr::lastTimeMarker() const
{
    return min
    (
        min
        (
            radialModel_->lastTimeMarker(), 
            azimuthalModel_->lastTimeMarker()
        ), 
        axialModel_->lastTimeMarker()
    );
}

// ************************************************************************* //

