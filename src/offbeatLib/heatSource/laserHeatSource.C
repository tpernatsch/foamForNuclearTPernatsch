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

#include "laserHeatSource.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchField.H"
#include "IndirectList.H"
#include "fvm.H"
#include "fvc.H"
    
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(laserHeatSource, 0);
    
    addToRunTimeSelectionTable
    (
        heatSource, 
        laserHeatSource, 
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::laserHeatSource::calcAddressing
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

Foam::laserHeatSource::laserHeatSource
(
    const fvMesh& mesh,
    const materials& materials,
    const dictionary& heatSourceDict,
    const bool timeDependent
)
:
    constantHeatSource(mesh, materials, heatSourceDict),
    I_
    (
        IOobject
        (
            "I",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    alphaAbs_(readScalar(heatSourceDict.lookup("alphaAbsorption"))),
    alphaAbsField_
    (
        IOobject
        (
            "alphaAbs",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless/dimLength, alphaAbs_),
        "zeroGradient"
    ),
#ifdef OPENFOAMFOUNDATION
    useEffectiveAlpha_(heatSourceDict.lookup<bool>("useEffectiveAlpha")),
    beamDirection_(heatSourceDict.lookup<vector>("beamDirection")),
#elif OPENFOAMESI
    useEffectiveAlpha_(heatSourceDict.get<bool>("useEffectiveAlpha")),
    beamDirection_(heatSourceDict.get<vector>("beamDirection")),
#endif
    beamDirectionField_
    (
        IOobject
        (
            "n",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("", dimless, beamDirection_),
        "zeroGradient"
    ),
    addr_()
{

    const globalOptions& globalOpt
    (mesh_.lookupObject<globalOptions>("globalOptions"));

    calcAddressing(heatSourceDict);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::laserHeatSource::~laserHeatSource()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::laserHeatSource::correct()
{
    const globalOptions& globalOpt
    (mesh_.lookupObject<globalOptions>("globalOptions"));

    if(mesh_.topoChanging())
    {
        calcAddressing(heatSourceDict_);
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
    )
    {
        if ( useEffectiveAlpha_ )
        {
            // Correct effective absorption XS due to porosity redistribution
            correctEffectiveAbsorption();
        }

        // Solve Equation for Laser Intensity attenuation
        surfaceScalarField flux = fvc::interpolate( beamDirectionField_ ) &  mesh_.Sf();

        fvScalarMatrix IEqn
        (
            fvm::div(flux, I_, "div(I)")
          + fvm::Sp(alphaAbsField_, I_)
        );

        IEqn.solve();

        // Build the resulting heat source
        scalarField& Qi = Q_.ref();

        forAll(addr_, addrI)
        {
            const label cellI = addr_[addrI];
            
            Qi[cellI] = I_[cellI] * alphaAbsField_[cellI];
        }

        // Reshape power profile due to porosity redistribution
        reshapePowerPorosity();
        
    }
    else
    {
        currentTime_ = mesh_.time().value();
    }

    Q_.correctBoundaryConditions();
}


Foam::scalar Foam::laserHeatSource::predictNextPowerDensity() const
{
    return averagePowerDensity();
}


Foam::scalar Foam::laserHeatSource::averagePowerDensity() const
{
    scalarField Vi(mesh_.V(), addr_);
    scalarField Qi(Q_, addr_);

    return gSum(Vi*Qi)/gSum(Vi);
}

void Foam::laserHeatSource::correctEffectiveAbsorption()
{
    if
    (
        mesh_.foundObject<volScalarField>("porosity")
    )
    {
        scalarField& alphaI = alphaAbsField_.ref();

        const scalarField& porosity = 
            mesh_.lookupObject<volScalarField>("porosity").internalField();

        // Reduced porosity field (only power producing cells)
        scalarField porosityFuelMat(porosity, addr_);

        // Reduced volumes field (only power producing cells)
        scalarField Vi(mesh_.V().field(), addr_);

        // Compute average porosity
        scalar avgPorosity = gSum(porosityFuelMat*Vi)/gSum(Vi);

        forAll(addr_, addrI)
        {
            const label cellI = addr_[addrI];

            // Get local value of porosity (bounded in [0,1])
            const scalar pI = max(min(porosity[cellI], 1.0), 0.0);

            // scalar PiLim = min(pI, 0.95); 
            scalar PiLim = min(pI, 1.0); 

            const scalar fp = (1-PiLim)/(1+0.5*PiLim);
            // const scalar fp = (1-PiLim);
            // const scalar fp = max(0.7-PiLim, 0.0);

            // Scale factor for volumetric power
            // alphaI[cellI] = alphaAbs_ * (1-pI)/(1-avgPorosity);
            alphaI[cellI] = alphaAbs_ * fp;
        }
    }
}

void Foam::laserHeatSource::reshapePowerPorosity()
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

        // Reference to porosity field
        const scalarField& porosity = 
            mesh_.lookupObject<volScalarField>("porosity").internalField();

        // Reduced porosity field (only power producing cells)
        scalarField porosityFuelMat(porosity, addr_);

        // Reduced volumes field (only power producing cells)
        scalarField Vi(mesh_.V().field(), addr_);

        // Compute average porosity
        avgPorosity = gSum(porosityFuelMat*Vi)/gSum(Vi);

        forAll(addr_, addrI)
        {
            const label cellI = addr_[addrI];

            // Get local value of porosity (bounded in [0,1])
            const scalar pI = max(min(porosity[cellI], 1.0), 0.0);

            // Scale factor for volumetric power
            Qi[cellI] *= (1-pI)/(1-avgPorosity);
        }
    }
}


// ************************************************************************* //
