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

#if defined __has_include
#  if __has_include(<commDataLayer.H>)
#    include <commDataLayer.H>
#    define isCommDataLayerIncluded
#  endif
#endif

#ifdef isCommDataLayerIncluded

#include "commDataLayer.H"
#include "fmiAxialProfileFvPatchScalarField.H"
#include "scalarFieldFieldINew.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

namespace Foam
{

//  * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

fmiAxialProfileFvPatchScalarField::
fmiAxialProfileFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedTemperatureFvPatchScalarField(p, iF),
    axialLocationsNameFromFMU_(),
    axialProfileNameFromFMU_(),
    zMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    zRef_(),
    dict_()
{}

fmiAxialProfileFvPatchScalarField::
fmiAxialProfileFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict,
    const bool valueRequired
)
:
    fixedTemperatureFvPatchScalarField(p, iF, dict, valueRequired),
    axialLocationsNameFromFMU_(dict.subDict("axialProfileDict").lookup("axialLocationsNameFromFMU")),
    axialProfileNameFromFMU_(dict.subDict("axialProfileDict").lookup("axialProfileNameFromFMU")),
    zMethod_(interpolateTableBase::interpolationMethodNames_[
        dict.subDict("axialProfileDict").lookupOrDefault<word>("axialInterpolationMethod", "linear")
    ]),
    zRef_(),
    dict_(dict.subDict("axialProfileDict"))
{
    scalarField axialLoc(0);
    scalarField profileData(dict.subDict("axialProfileDict").get<scalarField>("initialData"));

    // Set user defined or default axial profile
    if (dict.subDict("axialProfileDict").found("initialAxialLocations"))
    {
        axialLoc = dict.subDict("axialProfileDict").get<scalarField>("initialAxialLocations");
    }
    else
    {
        // forAll(profileData, locI)
        // {
            axialLoc.append(0/**(scalar)locI/((scalar)profileData.size()-1.0)*/);
        // }
    }

    // Stringify the intitial profile data from the user
    word initAxialLoc("");
    word initProfileData("");
    forAll(axialLoc, locI)
    {
        initAxialLoc += std::to_string(axialLoc[locI]) + " ";
        initProfileData += std::to_string(profileData[locI]) + " ";
    }

    // Communicating with the FMU
    commDataLayer& data = commDataLayer::New(this->db().time());

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

    // Face centers of the patch in the ref configuration
    const fvMesh& referenceMesh =
        patch().boundaryMesh().mesh().foundObject<fvMesh>("referenceMesh")
            ? patch().boundaryMesh().mesh().lookupObject<fvMesh>("referenceMesh")
            : patch().boundaryMesh().mesh();

    label patchID(patch().index());

    // Set main pin direction
    const globalOptions& globalOpt
    (
        patch().boundaryMesh().mesh().lookupObject<globalOptions>
        ("globalOptions")
    );

    vector pinDirection_ = (globalOpt.pinDirection());

    zRef_ = referenceMesh.boundaryMesh()[patchID].faceCentres() & pinDirection_;

    // Test if possible to build interpolation table
    if (axialLoc.size() >= 2)
    {
        // Build interpolation table
        scalarInterpolateTable zTable(axialLoc, profileData, zMethod_);

        // Loop over all the face centers of the boundary patch and update
        // the imposed patch temperature field
        forAll(zRef_, i)
        {
            T_[i] = zTable(zRef_[i]);
        }
    }
    else
    {
        forAll(zRef_, i)
        {
            T_[i] = profileData[0];
        }
    }

    fixedTemperatureFvPatchScalarField::updateSurfaceTemperature();

    fvPatchField<scalar>::updateCoeffs();
}


fmiAxialProfileFvPatchScalarField::
fmiAxialProfileFvPatchScalarField
(
    const fmiAxialProfileFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedTemperatureFvPatchScalarField(ptf, p, iF, mapper),
    axialLocationsNameFromFMU_(ptf.axialLocationsNameFromFMU_),
    axialProfileNameFromFMU_(ptf.axialProfileNameFromFMU_),
    zMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    zRef_(ptf.zRef_),
    dict_()
{}


fmiAxialProfileFvPatchScalarField::
fmiAxialProfileFvPatchScalarField
(
    const fmiAxialProfileFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedTemperatureFvPatchScalarField(ptf, iF),
    axialLocationsNameFromFMU_(ptf.axialLocationsNameFromFMU_),
    axialProfileNameFromFMU_(ptf.axialProfileNameFromFMU_),
    zMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    zRef_(ptf.zRef_),
    dict_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::fmiAxialProfileFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Communicating with the FMU
    commDataLayer& data = commDataLayer::New(this->db().time());

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
    scalarField axialProfile(0);
    std::stringstream ss(axialProfileFromFMU);
    while (getline(ss, value, ' '))
    {
        axialProfile.append(std::stod(value));
    }
    ss.str("");
    ss.clear();

    // Destringify axial locations
    scalarField axialLoc(0);
    ss << axialLocFromFMU;
    while (getline(ss, value, ' '))
    {
        axialLoc.append(std::stod(value));
    }

    if (axialLoc.size() >= 2)
    {
        // Build interpolation table
        scalarInterpolateTable zTable(axialLoc, axialProfile, zMethod_);

        // Loop over all the face centers of the boundary patch and update
        // the imposed patch temperature field
        forAll(zRef_, i)
        {
            T_[i] = zTable(zRef_[i]);
        }
    }
    else
    {
        forAll(zRef_, i)
        {
            T_[i] = axialProfile[0];
        }
    }

    fixedTemperatureFvPatchScalarField::updateSurfaceTemperature();

    fvPatchField<scalar>::updateCoeffs();

    dict_.set("initialData", axialProfile);
}

void Foam::fmiAxialProfileFvPatchScalarField::write(Ostream& os) const
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
    fmiAxialProfileFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif