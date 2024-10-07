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

#include "timeDependentAxialProfileFvPatchScalarField.H"
#include "scalarFieldFieldINew.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

namespace Foam
{
 
//  * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
 
timeDependentAxialProfileFvPatchScalarField::
timeDependentAxialProfileFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedTemperatureFvPatchScalarField(p, iF),
    zValues_(),
    timeValues_(),
    profileData_(),
    zMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    tMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    table_(timeValues_, profileData_, tMethod_),
    zRef_(),
    dict_()
{}

timeDependentAxialProfileFvPatchScalarField::
timeDependentAxialProfileFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict,
    const bool valueRequired
)
:
    fixedTemperatureFvPatchScalarField(p, iF, dict, valueRequired),
    zValues_(dict.subDict("axialProfileDict").lookup("axialLocations")),
    timeValues_(dict.subDict("axialProfileDict").lookup("timePoints")),
    profileData_
    (
        PtrList<scalarField>
        (
            dict.subDict("axialProfileDict").lookup("data"), 
            scalarFieldFieldINew()
        )
    ),
    zMethod_(interpolateTableBase::interpolationMethodNames_[
        dict.subDict("axialProfileDict").lookupOrDefault<word>("axialInterpolationMethod", "linear")
    ]),
    tMethod_(interpolateTableBase::interpolationMethodNames_[
        dict.subDict("axialProfileDict").lookupOrDefault<word>("timeInterpolationMethod", "linear")
    ]),
    table_(timeValues_, profileData_, tMethod_),
    zRef_(),
    dict_(dict.subDict("axialProfileDict"))
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

    // Set main pin direction
    vector pinDirection = db().foundObject<globalOptions>("globalOptions")?
    db().lookupObject<globalOptions>("globalOptions").pinDirection()
    : Foam::vector(0, 0, 1);

    zRef_ = referenceMesh.boundaryMesh()[patchID].faceCentres() & pinDirection;
 }
 
 
timeDependentAxialProfileFvPatchScalarField::
timeDependentAxialProfileFvPatchScalarField
(
    const timeDependentAxialProfileFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper,
    const bool mappingRequired
)
:
    fixedTemperatureFvPatchScalarField(ptf, p, iF, mapper, mappingRequired),
    zValues_(ptf.zValues_),
    timeValues_(ptf.timeValues_),
    profileData_(ptf.profileData_),
    zMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    tMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    table_(timeValues_, profileData_, tMethod_),
    zRef_(ptf.zRef_),
    dict_()
{}
 
 
timeDependentAxialProfileFvPatchScalarField::
timeDependentAxialProfileFvPatchScalarField
(
    const timeDependentAxialProfileFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedTemperatureFvPatchScalarField(ptf, iF),
    zValues_(ptf.zValues_),
    timeValues_(ptf.timeValues_),
    profileData_(ptf.profileData_),
    zMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    tMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    table_(timeValues_, profileData_, tMethod_),
    zRef_(ptf.zRef_),
    dict_()
{}
 
 
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
 

void Foam::timeDependentAxialProfileFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Value of current time
    const scalar &t = this ->db().time().value();

    // Create a scalar field for the profile values at time t
    scalarField zData
    (
        table_
        (
            patch().boundaryMesh().mesh().time().timeToUserTime(t)
        )
    );

    // Create interpolation table     
    scalarInterpolateTable zTable(zValues_, zData, zMethod_);

    // Loop over all the face centers of the boundary patch and update
    //  the imposed patch temperature field
    forAll(zRef_, i)
    {
        T_[i] = zTable(zRef_[i]);
    }

    fixedTemperatureFvPatchScalarField::updateSurfaceTemperature();

    fvPatchField<scalar>::updateCoeffs();
}
 
void Foam::timeDependentAxialProfileFvPatchScalarField::write(Ostream& os) const
{
    fixedTemperatureFvPatchScalarField::write(os);
#ifdef OPENFOAMFOUNDATION    
    writeEntry(os, "axialProfileDict", dict_);
#elif OPENFOAMESI
    os.writeEntry("axialProfileDict", dict_);
#endif    
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    timeDependentAxialProfileFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam