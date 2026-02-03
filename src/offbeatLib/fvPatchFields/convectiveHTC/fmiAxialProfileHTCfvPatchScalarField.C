/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#if defined __has_include
#  if __has_include(<commDataLayer.H>)
#    include <commDataLayer.H>
#    define isCommDataLayerIncluded
#  endif
#endif

#ifdef isCommDataLayerIncluded

#include "commDataLayer.H"
#include "fmiAxialProfileHTCfvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "scalarFieldFieldINew.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fmiAxialProfileHTCfvPatchScalarField::
fmiAxialProfileHTCfvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvPatchField<scalar>(p, iF),
    axialLocationsNameFromFMU_(),
    TaxialProfileNameFromFMU_(),
    HaxialProfileNameFromFMU_(),
    zMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    dict_(),
    kappaName_("k"),
    T0_(p.size(), 0.0),
    h_(p.size(), 0.0),
    alpha_(p.size(), 0.0)
{}


fmiAxialProfileHTCfvPatchScalarField::
fmiAxialProfileHTCfvPatchScalarField
 (
     const fvPatch& p,
     const DimensionedField<scalar, volMesh>& iF,
     const dictionary& dict,
     const bool valueRequired
 )
 :
    fvPatchField<scalar>(p, iF, dict, valueRequired),
    axialLocationsNameFromFMU_(dict.subDict("axialProfileDict").lookup("axialLocationsNameFromFMU")),
    TaxialProfileNameFromFMU_(dict.subDict("axialProfileDict").lookup("TaxialProfileNameFromFMU")),
    HaxialProfileNameFromFMU_(dict.subDict("axialProfileDict").lookup("HaxialProfileNameFromFMU")),
    zMethod_
    (
        interpolateTableBase::interpolationMethodNames_
        [
            dict.subDict("axialProfileDict").lookupOrDefault<word>
            ("axialInterpolationMethod", "linear")
        ]
    ),
    dict_(dict.subDict("axialProfileDict")),
    kappaName_(dict.lookupOrDefault<word>("kappa", "k")),
    T0_(p.size(), 0.0),
    h_(p.size(), 0.0),
    alpha_(p.size(), 0.0)
{
    scalarField axialLoc(0);
    scalarField TprofileData(dict.subDict("axialProfileDict").get<scalarField>("TinitialData"));
    scalarField HprofileData(dict.subDict("axialProfileDict").get<scalarField>("HinitialData"));

    // Set user defined or default axial profile
    if (dict.subDict("axialProfileDict").found("initialAxialLocations"))
    {
        axialLoc = dict.subDict("axialProfileDict").get<scalarField>("initialAxialLocations");
    }
    else
    {
        // forAll(TprofileData, locI)
        // {
            axialLoc.append(0/**(scalar)locI/((scalar)TprofileData.size()-1.0)*/);
        // }
    }

    // Stringify the intitial profile data from the user
    word initAxialLoc("");
    word initTProfileData("");
    word initHProfileData("");
    forAll(axialLoc, locI)
    {
        initAxialLoc += std::to_string(axialLoc[locI]) + " ";
        initTProfileData += std::to_string(TprofileData[locI]) + " ";
        initHProfileData += std::to_string(HprofileData[locI]) + " ";
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
        initTProfileData,
        TaxialProfileNameFromFMU_,
        commDataLayer::causality::in
    );
    data.storeObj
    (
        initHProfileData,
        HaxialProfileNameFromFMU_,
        commDataLayer::causality::in
    );
}


fmiAxialProfileHTCfvPatchScalarField::
fmiAxialProfileHTCfvPatchScalarField
(
    const fmiAxialProfileHTCfvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper,
    const bool mappingRequired
)
:
#ifdef OPENFOAMFOUNDATION
    fvPatchField<scalar>(ptf, p, iF, mapper, mappingRequired),
#elif OPENFOAMESI
    fvPatchField<scalar>(ptf, p, iF, mapper),
#endif
    axialLocationsNameFromFMU_(ptf.axialLocationsNameFromFMU_),
    TaxialProfileNameFromFMU_(ptf.TaxialProfileNameFromFMU_),
    HaxialProfileNameFromFMU_(ptf.HaxialProfileNameFromFMU_),
    zMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    dict_(),
    kappaName_(ptf.kappaName_),
    T0_(ptf.T0_),
    h_(ptf.h_),
    alpha_(ptf.alpha_)
{}


fmiAxialProfileHTCfvPatchScalarField::
fmiAxialProfileHTCfvPatchScalarField
(
    const fmiAxialProfileHTCfvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvPatchField<scalar>(tppsf, iF),
    axialLocationsNameFromFMU_(tppsf.axialLocationsNameFromFMU_),
    TaxialProfileNameFromFMU_(tppsf.TaxialProfileNameFromFMU_),
    HaxialProfileNameFromFMU_(tppsf.HaxialProfileNameFromFMU_),
    zMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    dict_(),
    kappaName_(tppsf.kappaName_),
    T0_(tppsf.T0_),
    h_(tppsf.h_),
    alpha_(tppsf.alpha_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void fmiAxialProfileHTCfvPatchScalarField::updateCoeffs()
{
    if (this->updated())
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
    const word HaxialProfileFromFMU = data.getObj<word>
    (
        HaxialProfileNameFromFMU_,
        commDataLayer::causality::in
    );

    // Destringify axial profile
    word value;
    scalarField HaxialProfile(0);
    std::stringstream ss(HaxialProfileFromFMU);
    while (getline(ss, value, ' '))
    {
        HaxialProfile.append(std::stod(value));
    }
    ss.str("");
    ss.clear();

    // Destringify axial profile
    scalarField axialLoc(0);
    ss << axialLocFromFMU;
    while (getline(ss, value, ' '))
    {
        axialLoc.append(std::stod(value));
    }

    scalarInterpolateTable h_zTable(axialLoc, HaxialProfile, zMethod_);

    const fvMesh& mesh = patch().boundaryMesh().mesh();
    const fvPatch& patch = this->patch();
    const scalarField& deltaCoeffs = patch.deltaCoeffs();
    scalarField alphaN(this->size(), 0.0);

    if (mesh.foundObject<volTensorField>(kappaName_))
    {
        const fvPatchField<tensor>& k = patch.lookupPatchField<volTensorField, tensor>
        (
            kappaName_
        );

        vectorField n = (patch.Sf() / patch.magSf());
        scalarField alphaP = ((n & k) & n)*deltaCoeffs;

        //- Loop over all the face centers
        forAll(this->patch().Cf(),faceI)
        {
            //- Z coordinate of the face i
            const scalar z = this->patch().Cf()[faceI].z();

            alphaN[faceI] = h_zTable(z);
        }

        alpha_ = alphaP / (alphaP + alphaN);
    }
    else
    {
        const fvPatchField<scalar>& k = patch.lookupPatchField<volScalarField, scalar>
        (
            kappaName_
        );

        scalarField alphaP = k*deltaCoeffs;

        //- Loop over all the face centers
        forAll(this->patch().Cf(),faceI)
        {
            //- Z coordinate of the face i
            const scalar z = this->patch().Cf()[faceI].z();

            alphaN[faceI] = h_zTable(z);
        }

        alpha_ = alphaP / (alphaP + alphaN);

        // Store h_ for write function
        h_ = alphaN;
    }

    fvPatchField<scalar>::updateCoeffs();
}


void fmiAxialProfileHTCfvPatchScalarField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
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
    const word TaxialProfileFromFMU = data.getObj<word>
    (
        TaxialProfileNameFromFMU_,
        commDataLayer::causality::in
    );

    // Destringify axial profile
    word value;
    scalarField TaxialProfile(0);
    std::stringstream ss(TaxialProfileFromFMU);
    while (getline(ss, value, ' '))
    {
        TaxialProfile.append(std::stod(value));
    }
    ss.str("");
    ss.clear();

    // Destringify axial profile
    scalarField axialLoc(0);
    ss << axialLocFromFMU;
    while (getline(ss, value, ' '))
    {
        axialLoc.append(std::stod(value));
    }

    scalarInterpolateTable T0_zTable(axialLoc, TaxialProfile, zMethod_);

    //- Loop over all the face centers
    forAll(this->patch().Cf(),faceI)
    {
        //- Z coordinate of the face i
        const scalar z = this->patch().Cf()[faceI].z();

        T0_[faceI] = T0_zTable(z);
    }

    fvPatchField<scalar>::operator=(
            alpha_*patchInternalField() + (1-alpha_)*T0_
    );

    fvPatchField<scalar>::evaluate();
}


tmp<Field<scalar> >
fmiAxialProfileHTCfvPatchScalarField::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<scalarField>(new scalarField(alpha_));
}


tmp<Field<scalar> >
fmiAxialProfileHTCfvPatchScalarField::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return (1-alpha_)*T0_;
}


tmp<Field<scalar> >
fmiAxialProfileHTCfvPatchScalarField::gradientInternalCoeffs() const
{
    const scalarField& deltaCoeffs = patch().deltaCoeffs();

    return (alpha_-1)*deltaCoeffs;
}


tmp<Field<scalar> >
fmiAxialProfileHTCfvPatchScalarField::gradientBoundaryCoeffs() const
{
    const scalarField& deltaCoeffs = patch().deltaCoeffs();

    return (1-alpha_)*deltaCoeffs*T0_;
}


void fmiAxialProfileHTCfvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
#ifdef OPENFOAMFOUNDATION

    if (!dict_.isNull())
    {
        writeEntry(os, "axialProfileDict", dict_);
    }

    writeEntryIfDifferent<word>(os, "kappa", "k", kappaName_);
    writeEntry(os, "T0"   , this->T0_);
    writeEntry(os, "h"    , this->h_);
    writeEntry(os, "value", *this);
#elif OPENFOAMESI

    if (!dict_.isNullDict())
    {
        os.writeEntry("axialProfileDict", dict_);
    }
    else
    {
        os.writeEntryIfDifferent<word>("kappa", "k", kappaName_);
        os.writeEntry("T0"   , T0_);
        os.writeEntry("h"    , h_);
        os.writeEntry("value", *this);
    }
#endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, fmiAxialProfileHTCfvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
