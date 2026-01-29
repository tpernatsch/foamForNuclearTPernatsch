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

#include "burnupDependentTabulatedRadialProfile.H"
#include "addToRunTimeSelectionTable.H"
#include "InterpolateTable.H"
#include "scalarFieldFieldINew.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(burnupDependentTabulatedRadialProfile, 0);

    addToRunTimeSelectionTable
    (
        radialProfile, 
        burnupDependentTabulatedRadialProfile, 
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::burnupDependentTabulatedRadialProfile::burnupDependentTabulatedRadialProfile
(
    const fvMesh& mesh,    
    const dictionary& dict,
    const labelList& addr,
    const vector& axialDirection,
    scalar rMax 
)
:
    radialProfile(mesh, dict, addr, axialDirection, rMax),
    rValues_(dict.lookup("radialLocations")),
    burnupValues_(dict.lookup("burnupPoints")),
    profileData_(PtrList<scalarField>(dict.lookup("data"), scalarFieldFieldINew())),
    rMethod_(interpolateTableBase::interpolationMethodNames_[
        dict.lookupOrDefault<word>("radialInterpolationMethod", "linear")
    ]),
    burnupMethod_(interpolateTableBase::interpolationMethodNames_[
        dict.lookupOrDefault<word>("burnupInterpolationMethod", "linear")
    ]),
    table_(burnupValues_, profileData_, burnupMethod_)
{
    if ((rValues_[0] != 0.0) || (rValues_[rValues_.size()-1] != 1.0))
    {
        FatalIOErrorInFunction(dict)
            << "Radial locations should be normalised from 0 to 1."
            << abort(FatalIOError);
    }
    
    //- Check that all profiles are normalised to 1
    forAll(profileData_, rowI)
    {
        scalarInterpolateTable rowTable
        (
            rValues_, 
            profileData_[rowI],
            rMethod_
        );
        
        scalar integral = rowTable.integral(2);
        
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

Foam::burnupDependentTabulatedRadialProfile::~burnupDependentTabulatedRadialProfile()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::burnupDependentTabulatedRadialProfile::correct()
{
    // Get average burnup
    const volScalarField& Bu = mesh_.lookupObject<volScalarField>("Bu");
    scalar avgBu;
    
    {
        scalarField Bui(IndirectList<scalar>(Bu.internalField(), addr_)());
        scalarField Vi(IndirectList<scalar>(mesh_.V(), addr_)());
        
        avgBu = gSum(Bui*Vi)/gSum(Vi);
    }
    
    if (debug)
    {
        Info<< "Radial burnup-dependent profile, average burnup is "
            << avgBu << " MWd/MT." << endl;
    }
    
    //- Interpolate burnup
    scalarField rData(table_(avgBu));
    
    scalarInterpolateTable rTable(rValues_, rData, rMethod_);
    
    scalarField r = this->radialLocations();
    
    forAll(addr_, i)
    {
        profile_[i] = rTable(r[i]);
    }
}

// ************************************************************************* //
