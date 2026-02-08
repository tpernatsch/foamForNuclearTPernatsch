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

#if defined __has_include
#  if __has_include(<commDataLayer.H>) 
#    include <commDataLayer.H>
#    define isCommDataLayerIncluded
#  endif
#endif

#ifdef isCommDataLayerIncluded

#include "commDataLayer.H"
#include "fmiTabulatedAxialProfile.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fmiTabulatedAxialProfile, 0);

    addToRunTimeSelectionTable
    (
        axialProfile, 
        fmiTabulatedAxialProfile, 
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fmiTabulatedAxialProfile::fmiTabulatedAxialProfile
(
    const fvMesh& mesh,    
    const dictionary& dict,
    const labelList& addr,
    const vector& axialDirection,
    scalar zMin,
    scalar zMax 
)
:
    axialProfile(mesh, dict, addr, axialDirection, zMin, zMax),
    axialLocationsNameFromFMU_(dict.lookup("axialLocationsNameFromFMU")),
    axialProfileNameFromFMU_(dict.lookup("axialProfileNameFromFMU")),
    zMethod_(interpolateTableBase::interpolationMethodNames_[
        dict.lookupOrDefault<word>("axialInterpolationMethod", "linear")
    ]),
    fmiState_
    (
        IOobject
        (
            "axialProfile",
            mesh_.time().timeName(),
            "uniform/fmiState",
            mesh_.time(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        )
    )
{
    scalarField axialLoc(0);
    scalarField profileData(0);

    const dictionary& subDict(fmiState_.subDict(axialProfileNameFromFMU_));

    if (subDict.found("values"))
    {
        profileData = subDict.get<scalarField>("values");
    }
    else
    {
        WarningIn("Foam::fmiTabulatedAxialProfile::fmiTabulatedAxialProfile")
            << nl
            << "    No 'values' has been found." << nl
            << "    Set normalized axial profile to 1.0." << nl
            << endl;
        profileData.append(1.0);
    }
    
    // Set user defined or default axial profile
    if (subDict.found("axialLocations"))
    {
        axialLoc = subDict.get<scalarField>("axialLocations");
    }
    else
    {
        WarningIn("Foam::fmiTabulatedAxialProfile::fmiTabulatedAxialProfile")
            << nl
            << "    No 'axialLocations' has been found." << nl
            << "    Set normalized axial profile from 0 to 1." << nl
            << endl;
        forAll(profileData, locI)
        {
            axialLoc.append((scalar)locI/((scalar)profileData.size()-1.0));
        }
    }

    // Stringify the initial profile data from the user
    word initAxialLoc("");
    word initProfileData("");
    forAll(axialLoc, locI)
    {
        initAxialLoc += std::to_string(axialLoc[locI]) + " ";
        initProfileData += std::to_string(profileData[locI]) + " ";
    }

    // Communicating with the FMU
    commDataLayer& data = commDataLayer::New(mesh_.time());

    // Store initial value
    data.storeObj
    (
        initAxialLoc,
        axialLocationsNameFromFMU_, 
        commDataLayer::causality::in
    );
    data.storeObj
    (
        initProfileData,
        axialProfileNameFromFMU_, 
        commDataLayer::causality::in
    );
    
    // Test if possible to build interpolation table
    if (axialLoc.size() >= 2)
    {
        scalarInterpolateTable zTable(axialLoc, profileData, zMethod_);
        
        if ((axialLoc[0] != 0.0) || (axialLoc[axialLoc.size()-1] != 1.0))
        {
            FatalIOErrorInFunction(fmiState_)
                << "Axial locations should be normalised from 0 to 1."
                << abort(FatalIOError);
        }
        
        // Check that profile is normalised to 1
        scalar integral = zTable.integral(1);
        
        if (mag(integral - 1) > 1e-4)
        {
            FatalIOErrorInFunction(fmiState_)
                << "Supplied profile is not normalised to 1. "
                << "Current integral is " << integral
                << endl << abort(FatalIOError);
        }
        
        // Check that profile contains no negative values
        if (max(profileData) < 0)
        {
            FatalIOErrorInFunction(fmiState_)
                << "Supplied profile contains negative values."
                << abort(FatalIOError);
        }
        
        //- Update profile
        scalarField z = this->axialLocations();
        
        forAll(addr_, i)
        {
            // Unnormalize and renormalize to the axialLocation from FMI
            // If zMin_ == zMin && zMax_ == zMax then there no change
            profile_[i] = zTable((z[i] * (zMax_ - zMin_) + zMin_ - zMin) / (zMax - zMin));
        }
    }
    else
    {
        forAll(addr_, i)
        {
            profile_[i] = profileData[0];
        }
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fmiTabulatedAxialProfile::~fmiTabulatedAxialProfile()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fmiTabulatedAxialProfile::correct()
{
    // Communicating with the FMU
    commDataLayer& data = commDataLayer::New(mesh_.time());

    // Get axial locations from the FMU
    const word axialLocFromFMU = data.getObj<word>
    (
        axialLocationsNameFromFMU_, 
        commDataLayer::causality::in
    );

    // Get axial profile from the FMU
    const word axialProfileFromFMU = data.getObj<word>
    (
        axialProfileNameFromFMU_, 
        commDataLayer::causality::in
    );
    
    // Destringify axial profile
    word value;
    scalarField profileData(0);
    std::stringstream ss(axialProfileFromFMU);
    while (getline(ss, value, ' '))
    {
        profileData.append(std::stod(value));
    }
    ss.str("");
    ss.clear();

    // Destringify and normalize axial locations
    // Axial locations are in absolute, not relative/normalized
    scalarField axialLoc(0);
    ss << axialLocFromFMU;
    while (getline(ss, value, ' '))
    {
        // axialLoc.append((std::stod(value) - zMin_) / (zMax_ - zMin_));
        axialLoc.append(std::stod(value));
    }
    const scalar zMin(axialLoc[0]);
    const scalar zMax(axialLoc[axialLoc.size()-1]);
    forAll(axialLoc, i)
    {
        axialLoc[i] = (axialLoc[i] - zMin) / (zMax - zMin);
    }

    // Test if possible to build interpolation table
    if (axialLoc.size() >= 2)
    {
        // Build interpolation table
        scalarInterpolateTable zTable(axialLoc, profileData, zMethod_);

        // Update profile in mesh
        scalarField z = this->axialLocations();
        forAll(addr_, i)
        {
            // Unnormalize and renormalize to the axialLocation from FMI
            // If zMin_ == zMin && zMax_ == zMax then there no change
            profile_[i] = zTable((z[i] * (zMax_ - zMin_) + zMin_ - zMin) / (zMax - zMin));
        }
    }
    else
    {
        forAll(addr_, i)
        {
            profile_[i] = profileData[0];
        }
    }


    //- Write FMI state
    if (mesh_.time().writeTime())
    {
        fmiState_.subDict(axialProfileNameFromFMU_).set("axialLocations", axialLoc);
        fmiState_.subDict(axialProfileNameFromFMU_).set("values", profileData);
    }
}

#endif

// ************************************************************************* //
