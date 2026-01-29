/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
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

#include "timeDependentAzimuthalProfileFvPatchScalarField.H"
#include "scalarFieldFieldINew.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

namespace Foam
{
 
//  * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
 
timeDependentAzimuthalProfileFvPatchScalarField::
timeDependentAzimuthalProfileFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedTemperatureFvPatchScalarField(p, iF),
    thetaValues_(),
    timeValues_(),
    profileData_(),
    thetaMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    tMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    table_(timeValues_, profileData_, tMethod_),
    thetaRef_(p.size(), 0.0),
    dict_()
{}

timeDependentAzimuthalProfileFvPatchScalarField::
timeDependentAzimuthalProfileFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict,
    const bool valueRequired
)
:
    fixedTemperatureFvPatchScalarField(p, iF, dict, valueRequired),
    thetaValues_(dict.subDict("azimuthalProfileDict").lookup("azimuthalLocations")),
    timeValues_(dict.subDict("azimuthalProfileDict").lookup("timePoints")),
    profileData_
    (
        PtrList<scalarField>
        (
            dict.subDict("azimuthalProfileDict").lookup("data"), 
            scalarFieldFieldINew()
        )
    ),
    thetaMethod_(interpolateTableBase::interpolationMethodNames_[
        dict.subDict("azimuthalProfileDict").lookupOrDefault<word>("azimuthalInterpolationMethod", "linear")
    ]),
    tMethod_(interpolateTableBase::interpolationMethodNames_[
        dict.subDict("azimuthalProfileDict").lookupOrDefault<word>("timeInterpolationMethod", "linear")
    ]),
    table_(timeValues_, profileData_, tMethod_),
    thetaRef_(p.size(), 0.0),
    dict_(dict.subDict("azimuthalProfileDict"))
{
    // Face centers of the patch in the ref configuration
    const fvMesh& referenceMesh = 
    ( 
    patch().boundaryMesh().mesh().foundObject<fvMesh>("referenceMesh")
    ) ?
    patch().boundaryMesh().mesh().lookupObject<fvMesh>("referenceMesh")
    :
    patch().boundaryMesh().mesh();

    label patchID(patch().index());

    // Compute theta of the patch faces
    vectorField faceCentres =  referenceMesh.boundaryMesh()[patchID].faceCentres();

    // Assumed pin direction along z
    forAll(faceCentres, i)
    {
        vector center(faceCentres[i]);

        // Coordinates of this face
        scalar x = center.x();
        scalar y = center.y();

        // Calculate the angle theta in degrees 
        thetaRef_[i] = atan2(y, x)* (180.0 / Foam::constant::mathematical::pi);
    }
 }
 
 
timeDependentAzimuthalProfileFvPatchScalarField::
timeDependentAzimuthalProfileFvPatchScalarField
(
    const timeDependentAzimuthalProfileFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper,
    const bool mappingRequired
)
:
    fixedTemperatureFvPatchScalarField(ptf, p, iF, mapper, mappingRequired),
    thetaValues_(ptf.thetaValues_),
    timeValues_(ptf.timeValues_),
    profileData_(ptf.profileData_),
    thetaMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    tMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    table_(timeValues_, profileData_, tMethod_),
    thetaRef_(ptf.thetaRef_),
    dict_()
{}
 
 
timeDependentAzimuthalProfileFvPatchScalarField::
timeDependentAzimuthalProfileFvPatchScalarField
(
    const timeDependentAzimuthalProfileFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedTemperatureFvPatchScalarField(ptf, iF),
    thetaValues_(ptf.thetaValues_),
    timeValues_(ptf.timeValues_),
    profileData_(ptf.profileData_),
    thetaMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    tMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    table_(timeValues_, profileData_, tMethod_),
    thetaRef_(ptf.thetaRef_),
    dict_()
{}
 
 
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
 

void Foam::timeDependentAzimuthalProfileFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Value of current time
    const scalar &t = this ->db().time().value();

    // Create a scalar field for the profile values at time t
    scalarField thetaData
    (
        table_
        (
            patch().boundaryMesh().mesh().time().timeToUserTime(t)
        )
    );

    // Create interpolation table     
    scalarInterpolateTable thetaTable(thetaValues_, thetaData, thetaMethod_);

    // Loop over all the face centers of the boundary patch and update
    //  the imposed patch temperature field
    forAll(thetaRef_, i)
    {
        T_[i] = thetaTable(thetaRef_[i]);
    }

    fixedTemperatureFvPatchScalarField::updateSurfaceTemperature();

    fvPatchField<scalar>::updateCoeffs();
}
 
void Foam::timeDependentAzimuthalProfileFvPatchScalarField::write(Ostream& os) const
{
    fixedTemperatureFvPatchScalarField::write(os);
#ifdef OPENFOAMFOUNDATION    
    writeEntry(os, "azimuthalProfileDict", dict_);
#elif OPENFOAMESI
    os.writeEntry("azimuthalProfileDict", dict_);
#endif    
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    timeDependentAzimuthalProfileFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam