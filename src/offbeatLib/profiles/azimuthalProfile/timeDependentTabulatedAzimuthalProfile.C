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

#include "timeDependentTabulatedAzimuthalProfile.H"
#include "addToRunTimeSelectionTable.H"
#include "InterpolateTable.H"
#include "scalarFieldFieldINew.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(timeDependentTabulatedAzimuthalProfile, 0);

    addToRunTimeSelectionTable
    (
        azimuthalProfile, 
        timeDependentTabulatedAzimuthalProfile, 
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeDependentTabulatedAzimuthalProfile::timeDependentTabulatedAzimuthalProfile
(
    const fvMesh& mesh,    
    const dictionary& dict,
    const labelList& addr,
    const vector& axialDirection,
    const vector& radialDirection
)
:
    azimuthalProfile(mesh, dict, addr, axialDirection, radialDirection),
    thetaValues_(dict.lookup("azimuthalLocations")),
    timeValues_(dict.lookup("timePoints")),
    profileData_(PtrList<scalarField>(dict.lookup("data"), scalarFieldFieldINew())),
    thetaMethod_(interpolateTableBase::interpolationMethodNames_[
        dict.lookupOrDefault<word>("azimuthalInterpolationMethod", "linear")
    ]),
    tMethod_(interpolateTableBase::interpolationMethodNames_[
        dict.lookupOrDefault<word>("timeInterpolationMethod", "linear")
    ]),
    table_(timeValues_, profileData_, tMethod_)
{
    if ((thetaValues_[0] != 0.0) || (thetaValues_[thetaValues_.size()-1] != 1.0))
    {
        FatalIOErrorInFunction(dict)
            << "Azimuthal locations should be normalised from 0 to 1."
            << abort(FatalIOError);
    }
    
    //- Check that all profiles are normalised to 1
    forAll(profileData_, rowI)
    {
        scalarInterpolateTable rowTable
        (
            thetaValues_, 
            profileData_[rowI],
            thetaMethod_
        );
        
        scalar integral = rowTable.integral(1);
        
        if (mag(integral - 1) > 1e-4)
        {
            FatalIOErrorInFunction(dict)
                << "Supplied profile is not normalised to 1 for table row "
                << rowI << ". Current integral is " << integral
                << endl << abort(FatalIOError);
        }
    }

    if (max(profileData_) < 0)
    {
        FatalIOErrorInFunction(dict)
            << "Supplied profile data contains negative values."
            << abort(FatalIOError);
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::timeDependentTabulatedAzimuthalProfile::~timeDependentTabulatedAzimuthalProfile()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::timeDependentTabulatedAzimuthalProfile::correct()
{
    //- Interpolate time
    scalar t = mesh_.time().value();
    scalarField thetaData(table_(mesh_.time().timeToUserTime(t)));
    
    scalarInterpolateTable thetaTable(thetaValues_, thetaData, thetaMethod_);
    
    //- TODO is it necessary to call azimuthalLocations every time?
    scalarField theta = this->azimuthalLocations();
    
    forAll(addr_, i)
    {
        profile_[i] = thetaTable(theta[i]);
    }
}

// ************************************************************************* //
