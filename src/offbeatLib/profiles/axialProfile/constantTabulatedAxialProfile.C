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

#include "constantTabulatedAxialProfile.H"
#include "addToRunTimeSelectionTable.H"
#include "InterpolateTables.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(constantTabulatedAxialProfile, 0);

    addToRunTimeSelectionTable
    (
        axialProfile, 
        constantTabulatedAxialProfile, 
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constantTabulatedAxialProfile::constantTabulatedAxialProfile
(
    const fvMesh& mesh,    
    const dictionary& dict,
    const labelList& addr,
    const vector& axialDirection,
    scalar zMin,
    scalar zMax 
)
:
    axialProfile(mesh, dict, addr, axialDirection, zMin, zMax)
{
    scalarField zValues(dict.lookup("axialLocations"));

    scalarField scalingFactorData(dict.found("axialScalingFactors") ?
        dict.lookup("axialScalingFactors") :
        scalarField(zValues.size(), 1.0));
    scalarField profileData(dict.lookup("data"));
    interpolateTableBase::interpolationMethod method =
        interpolateTableBase::interpolationMethodNames_[
            word(dict.lookup("axialInterpolationMethod"))];
    
    scalarInterpolateTable zTable(zValues, profileData, method);
    scalarInterpolateTable scalingFactorTable(
        zValues, scalingFactorData, method);
    
    if ((zValues[0] != 0.0) || (zValues[zValues.size()-1] != 1.0))
    {
        FatalIOErrorInFunction(dict)
            << "Axial locations should be normalised from 0 to 1."
            << abort(FatalIOError);
    }
    
    // Check that profile is normalised to 1
    scalar integral = zTable.integral(1);
    
    if (mag(integral - 1) > 1e-4)
    {
        FatalIOErrorInFunction(dict)
            << "Supplied profile is not normalised to 1. "
            << "Current integral is " << integral
            << endl << abort(FatalIOError);
    }
    
    // Check that profile contains no negative values
    if (max(profileData) < 0)
    {
        FatalIOErrorInFunction(dict)
            << "Supplied profile contains negative values."
            << abort(FatalIOError);
    }

    if (max(scalingFactorData) < 0)
    {
        FatalIOErrorInFunction(dict)
            << "Supplied axial scaling factor data contains negative values."
            << abort(FatalIOError);
    }
    
    //- Update profile
    scalarField z = this->axialLocations();
    
    forAll(addr_, i)
    {
        profile_[i] = zTable(z[i]);
        scalingFactorProfile_[i] = scalingFactorTable(z[i]);
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::constantTabulatedAxialProfile::~constantTabulatedAxialProfile()
{}

// ************************************************************************* //
