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

#include "constantTabulatedRadialProfile.H"
#include "addToRunTimeSelectionTable.H"
#include "InterpolateTables.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(constantTabulatedRadialProfile, 0);

    addToRunTimeSelectionTable
    (
        radialProfile, 
        constantTabulatedRadialProfile, 
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constantTabulatedRadialProfile::constantTabulatedRadialProfile
(
    const fvMesh& mesh,    
    const dictionary& dict,
    const labelList& addr,
    const vector& axialDirection,
    scalar rMax 
)
:
    radialProfile(mesh, dict, addr, axialDirection, rMax)
{
    scalarField rData(dict.lookup("radialLocations"));
    scalarField profileData(dict.lookup("data"));
    interpolateTableBase::interpolationMethod method =
        interpolateTableBase::interpolationMethodNames_[
            word(dict.lookup("radialInterpolationMethod"))];
    
    scalarInterpolateTable rTable(rData, profileData, method);
    
    if ((rData[0] != 0.0) || (rData[rData.size()-1] != 1.0))
    {
        FatalIOErrorInFunction(dict)
            << "Axial locations should be normalised from 0 to 1."
            << abort(FatalIOError);
    }
    
    // Check that profile is normalised to 1
    scalar integral = rTable.integral(2);
    
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
    scalarField r = this->radialLocations();
    
    forAll(addr_, i)
    {
        profile_[i] = rTable(r[i]);
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::constantTabulatedRadialProfile::~constantTabulatedRadialProfile()
{}

// ************************************************************************* //
