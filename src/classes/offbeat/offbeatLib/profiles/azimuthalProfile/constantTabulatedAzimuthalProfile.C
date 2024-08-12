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

#include "constantTabulatedAzimuthalProfile.H"
#include "addToRunTimeSelectionTable.H"
#include "InterpolateTables.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(constantTabulatedAzimuthalProfile, 0);

    addToRunTimeSelectionTable
    (
        azimuthalProfile, 
        constantTabulatedAzimuthalProfile, 
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constantTabulatedAzimuthalProfile::constantTabulatedAzimuthalProfile
(
    const fvMesh& mesh,    
    const dictionary& dict,
    const labelList& addr,
    const vector& axialDirection,
    const vector& radialDirection
)
:
    azimuthalProfile(mesh, dict, addr, axialDirection, radialDirection)
{
    scalarField thetaData(dict.lookup("azimuthalLocations"));
    scalarField profileData(dict.lookup("data"));
    interpolateTableBase::interpolationMethod method =
        interpolateTableBase::interpolationMethodNames_[
            word(dict.lookup("azimuthalInterpolationMethod"))];
    
    scalarInterpolateTable thetaTable(thetaData, profileData, method);
    
    if ((thetaData[0] != 0.0) || (thetaData[thetaData.size()-1] != 1.0))
    {
        FatalIOErrorInFunction(dict)
            << "Axial locations should be normalised from 0 to 1."
            << abort(FatalIOError);
    }
    
    // Check that profile is normalised to 1
    scalar integral = thetaTable.integral(1);
    
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
    
    //- Update profile
    scalarField theta = this->azimuthalLocations();
    
    forAll(addr_, i)
    {
        profile_[i] = thetaTable(theta[i]);
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::constantTabulatedAzimuthalProfile::~constantTabulatedAzimuthalProfile()
{}

// ************************************************************************* //
