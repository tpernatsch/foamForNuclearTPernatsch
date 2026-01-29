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

#include "timeDependentAxialProfileHTCfvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "scalarFieldFieldINew.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

timeDependentAxialProfileHTCfvPatchScalarField::
timeDependentAxialProfileHTCfvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvPatchField<scalar>(p, iF),
    zValues_(),
    timeValues_(),
    hProfileData_(),
    T0ProfileData_(),
    zMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    tMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    hTable_(timeValues_, hProfileData_, tMethod_),
    T0Table_(timeValues_, T0ProfileData_, tMethod_),
    dict_(),
    kappaName_("k"),
    T0_(p.size(), 0.0),
    h_(p.size(), 0.0),
    alpha_(p.size(), 0.0)
{}


timeDependentAxialProfileHTCfvPatchScalarField::
timeDependentAxialProfileHTCfvPatchScalarField
 (
     const fvPatch& p,
     const DimensionedField<scalar, volMesh>& iF,
     const dictionary& dict,
     const bool valueRequired
 )
 :
    fvPatchField<scalar>(p, iF, dict, valueRequired),
    zValues_(dict.subDict("axialProfileDict").lookup("axialLocations")),
    timeValues_(dict.subDict("axialProfileDict").lookup("timePoints")),
    hProfileData_
    (
        PtrList<scalarField>(dict.subDict("axialProfileDict").lookup("hData"), 
            scalarFieldFieldINew())
    ),
    T0ProfileData_
    (
        PtrList<scalarField>(dict.subDict("axialProfileDict").lookup("T0Data"), 
            scalarFieldFieldINew())
    ),
    zMethod_
    (
        interpolateTableBase::interpolationMethodNames_
        [
            dict.subDict("axialProfileDict").lookupOrDefault<word>
            ("axialInterpolationMethod", "linear")
        ]
    ),
    tMethod_
    (
        interpolateTableBase::interpolationMethodNames_
        [
            dict.subDict("axialProfileDict").lookupOrDefault<word>
            ("timeInterpolationMethod", "linear")
        ]
    ),
    hTable_(timeValues_, hProfileData_, tMethod_),
    T0Table_(timeValues_, T0ProfileData_, tMethod_),
    dict_(dict.subDict("axialProfileDict")),
    kappaName_(dict.lookupOrDefault<word>("kappa", "k")),
    T0_(p.size(), 0.0),
    h_(p.size(), 0.0),
    alpha_(p.size(), 0.0)
{}

timeDependentAxialProfileHTCfvPatchScalarField::
timeDependentAxialProfileHTCfvPatchScalarField
(
    const timeDependentAxialProfileHTCfvPatchScalarField& ptf,
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
    zValues_(ptf.zValues_),
    timeValues_(ptf.timeValues_),
    hProfileData_(ptf.hProfileData_),
    T0ProfileData_(ptf.T0ProfileData_),
    zMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    tMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    hTable_(timeValues_, hProfileData_, tMethod_),
    T0Table_(timeValues_, T0ProfileData_, tMethod_),
    dict_(),
    kappaName_(ptf.kappaName_),
    T0_(ptf.T0_),
    h_(ptf.h_),
    alpha_(ptf.alpha_)
{}


timeDependentAxialProfileHTCfvPatchScalarField::
timeDependentAxialProfileHTCfvPatchScalarField
(
    const timeDependentAxialProfileHTCfvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvPatchField<scalar>(tppsf, iF),
    zValues_(tppsf.zValues_),
    timeValues_(tppsf.timeValues_),
    hProfileData_(tppsf.hProfileData_),
    T0ProfileData_(tppsf.T0ProfileData_),
    zMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    tMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    hTable_(timeValues_, hProfileData_, tMethod_),
    T0Table_(timeValues_, T0ProfileData_, tMethod_),
    dict_(),
    kappaName_(tppsf.kappaName_),
    T0_(tppsf.T0_),
    h_(tppsf.h_),
    alpha_(tppsf.alpha_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void timeDependentAxialProfileHTCfvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Current time 
    const scalar &t = this ->db().time().value();
    
    // Create interpolate table for h 
    scalarField h_zData
    (
        hTable_
        (
            patch().boundaryMesh().mesh().time().timeToUserTime(t)
        )
    );

    scalarInterpolateTable h_zTable(zValues_, h_zData, zMethod_);

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


void timeDependentAxialProfileHTCfvPatchScalarField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    // Current time 
    const scalar &t = this ->db().time().value();

    // Create interpolate table for T0_
    scalarField T0_zData
    (
        T0Table_
        (
            patch().boundaryMesh().mesh().time().timeToUserTime(t)
        )
    );
    
    scalarInterpolateTable T0_zTable(zValues_, T0_zData, zMethod_);

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
timeDependentAxialProfileHTCfvPatchScalarField::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<scalarField>(new scalarField(alpha_));
}


tmp<Field<scalar> >
timeDependentAxialProfileHTCfvPatchScalarField::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return (1-alpha_)*T0_;
}


tmp<Field<scalar> >
timeDependentAxialProfileHTCfvPatchScalarField::gradientInternalCoeffs() const
{
    const scalarField& deltaCoeffs = patch().deltaCoeffs();

    return (alpha_-1)*deltaCoeffs;
}


tmp<Field<scalar> >
timeDependentAxialProfileHTCfvPatchScalarField::gradientBoundaryCoeffs() const
{
    const scalarField& deltaCoeffs = patch().deltaCoeffs();

    return (1-alpha_)*deltaCoeffs*T0_;
}


void timeDependentAxialProfileHTCfvPatchScalarField::write(Ostream& os) const
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
    
    os.writeEntryIfDifferent<word>("kappa", "k", kappaName_);
    T0_.writeEntry("T0"   , os);
    h_.writeEntry("h"    , os);
    this->writeEntry("value", os);
#endif  
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, timeDependentAxialProfileHTCfvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
