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

#include "xzTemperatureProfileFvPatchScalarField.H"
#include "scalarFieldFieldINew.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

namespace Foam
{
 
//  * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
 
xzTemperatureProfileFvPatchScalarField::
xzTemperatureProfileFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedTemperatureFvPatchScalarField(p, iF),
    zValues_(),
    xValues_(),
    profileData_(),
    zMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    xMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    table_(xValues_, profileData_, xMethod_),
    dict_()
{}

xzTemperatureProfileFvPatchScalarField::
xzTemperatureProfileFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict,
    const bool valueRequired
)
:
    fixedTemperatureFvPatchScalarField(p, iF, dict, valueRequired),
    zValues_(dict.subDict("axialProfileDict").lookup("zLocations")),
    xValues_(dict.subDict("axialProfileDict").lookup("xLocations")),
    profileData_
    (
        PtrList<scalarField>
        (
            dict.subDict("axialProfileDict").lookup("data"), 
            scalarFieldFieldINew()
        )
    ),
    zMethod_(interpolateTableBase::interpolationMethodNames_[
        dict.subDict("axialProfileDict").lookupOrDefault<word>("zInterpolationMethod", "linear")
    ]),
    xMethod_(interpolateTableBase::interpolationMethodNames_[
        dict.subDict("axialProfileDict").lookupOrDefault<word>("xInterpolationMethod", "linear")
    ]),
    table_(xValues_, profileData_, xMethod_),
    dict_(dict.subDict("axialProfileDict"))
{}
 
 
xzTemperatureProfileFvPatchScalarField::
xzTemperatureProfileFvPatchScalarField
(
    const xzTemperatureProfileFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedTemperatureFvPatchScalarField(ptf, p, iF, mapper),
    zValues_(ptf.zValues_),
    xValues_(ptf.xValues_),
    profileData_(ptf.profileData_),
    zMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    xMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    table_(xValues_, profileData_, xMethod_),
    dict_()
{}
 
 
xzTemperatureProfileFvPatchScalarField::
xzTemperatureProfileFvPatchScalarField
(
    const xzTemperatureProfileFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedTemperatureFvPatchScalarField(ptf, iF),
    zValues_(ptf.zValues_),
    xValues_(ptf.xValues_),
    profileData_(ptf.profileData_),
    zMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    xMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    table_(xValues_, profileData_, xMethod_),
    dict_()
{}
 
 
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
 

void Foam::xzTemperatureProfileFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Take reference to mesh
    const fvMesh& referenceMesh = 
    ( 
    patch().boundaryMesh().mesh().foundObject<fvMesh>("referenceMesh")
    ) ?
    patch().boundaryMesh().mesh().lookupObject<fvMesh>("referenceMesh")
    :
    patch().boundaryMesh().mesh();

    // Label of current patch
    label patchID(patch().index());

    // Save face center vector field
    vectorField faceCentres(
        referenceMesh.boundaryMesh()[patchID].faceCentres());

    // For all the elements of the patch
    forAll(faceCentres, i)
    {
        // Take the x and z coordinates of the element
        scalar x(faceCentres[i].x());
        scalar z(faceCentres[i].z());

        // Create a scalar field for the profile values at location x
        scalarField zData(table_(x));

        // Create a z-dependent table at location x
        scalarInterpolateTable zTable(zValues_, zData, zMethod_);

        T_[i] = zTable(z);
    }  

    fixedTemperatureFvPatchScalarField::updateSurfaceTemperature();

    fvPatchField<scalar>::updateCoeffs();
}
 
void Foam::xzTemperatureProfileFvPatchScalarField::write(Ostream& os) const
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
    xzTemperatureProfileFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam