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

#include "azimuthalProfileFvPatchScalarField.H"
#include "scalarFieldFieldINew.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

namespace Foam
{
 
//  * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
 
azimuthalProfileFvPatchScalarField::
azimuthalProfileFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedTemperatureFvPatchScalarField(p, iF),
    thetaValues_(),
    profileData_(p.size(), 290),
    thetaMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    table_(thetaValues_, profileData_, thetaMethod_),
    thetaRef_(p.size(), 0),
    dict_()
{}

azimuthalProfileFvPatchScalarField::
azimuthalProfileFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict,
    const bool valueRequired
)
:
    fixedTemperatureFvPatchScalarField(p, iF, dict, valueRequired),
    thetaValues_(dict.subDict("azimuthalProfileDict").lookup("azimuthalLocations")),
    profileData_(dict.subDict("azimuthalProfileDict").lookup("data")),
    thetaMethod_(interpolateTableBase::interpolationMethodNames_[
        dict.subDict("azimuthalProfileDict").lookupOrDefault<word>("azimuthalInterpolationMethod", "linear")
    ]),
    table_(thetaValues_, profileData_, thetaMethod_),
    thetaRef_(p.size(), 0),
    dict_(dict.subDict("azimuthalProfileDict"))
{
    //- Face centers of the patch in the ref configuration
    const fvMesh& referenceMesh = 
    ( 
    patch().boundaryMesh().mesh().foundObject<fvMesh>("referenceMesh")
    ) ?
    patch().boundaryMesh().mesh().lookupObject<fvMesh>("referenceMesh")
    :
    patch().boundaryMesh().mesh();

    label patchID(patch().index());

//- Set main pin direction
    const globalOptions& globalOpt
    (
        patch().boundaryMesh().mesh().lookupObject<globalOptions>
        ("globalOptions")
    );

    vector pinDirection_ = (globalOpt.pinDirection());

    // Compute theta of the patch faces
    vectorField faceCentres =  referenceMesh.boundaryMesh()[patchID].faceCentres();

    // Assumed pin direction along z
    forAll(faceCentres, i)
    {
        vector center(faceCentres[i]);

        // Coordinates of this face
        scalar x = center.x();
        scalar y = center.y();

        // Calculate the angle theta
        thetaRef_[i] = atan2(y, x)* (180.0 / Foam::constant::mathematical::pi);
    }

    //- Loop over all the face centers of the boundary patch and update
    //  the imposed patch temperature field
    forAll(thetaRef_, i)
    {
        T_[i] = table_(thetaRef_[i]);
    }

    fixedTemperatureFvPatchScalarField::updateSurfaceTemperature();

    fvPatchField<scalar>::updateCoeffs();
 }
 
 
azimuthalProfileFvPatchScalarField::
azimuthalProfileFvPatchScalarField
(
    const azimuthalProfileFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper,
    const bool mappingRequired
)
:
    fixedTemperatureFvPatchScalarField(ptf, p, iF, mapper, mappingRequired),
    thetaValues_(ptf.thetaValues_),
    profileData_(ptf.profileData_),
    thetaMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    table_(thetaValues_, profileData_, thetaMethod_),
    thetaRef_(ptf.thetaRef_),
    dict_()
{}
 
 
azimuthalProfileFvPatchScalarField::
azimuthalProfileFvPatchScalarField
(
    const azimuthalProfileFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedTemperatureFvPatchScalarField(ptf, iF),
    thetaValues_(ptf.thetaValues_),
    profileData_(ptf.profileData_),
    thetaMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    table_(thetaValues_, profileData_, thetaMethod_),
    thetaRef_(ptf.thetaRef_),
    dict_()
{}
 
 
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
 

void Foam::azimuthalProfileFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    fvPatchField<scalar>::updateCoeffs();
}
 
void Foam::azimuthalProfileFvPatchScalarField::write(Ostream& os) const
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
    azimuthalProfileFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam