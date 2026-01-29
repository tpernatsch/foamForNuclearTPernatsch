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

#include "fixedHeatFluxFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fixedHeatFluxFvPatchScalarField::fixedHeatFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    q_(p.size(), 0.0),
    kappaName_("k")
{}


fixedHeatFluxFvPatchScalarField::fixedHeatFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF, dict),
    q_("heatFlux", dict, p.size()),
    kappaName_(dict.lookupOrDefault<word>("kappa", "k"))
{}


fixedHeatFluxFvPatchScalarField::fixedHeatFluxFvPatchScalarField
(
    const fixedHeatFluxFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    q_(mapper(ptf.q_)),
    kappaName_(ptf.kappaName_)
{}


fixedHeatFluxFvPatchScalarField::fixedHeatFluxFvPatchScalarField
(
    const fixedHeatFluxFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(tppsf, iF),
    q_(tppsf.q_),
    kappaName_(tppsf.kappaName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void fixedHeatFluxFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchScalarField::autoMap(m);

#ifdef OPENFOAMFOUNDATION
    m(q_, q_);
#elif OPENFOAMESI
    q_.autoMap(m);
#endif
}


void fixedHeatFluxFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchScalarField::rmap(ptf, addr);

    const fixedHeatFluxFvPatchScalarField& dmptf =
        refCast<const fixedHeatFluxFvPatchScalarField>(ptf);

    q_.rmap(dmptf.q_, addr);
}


void fixedHeatFluxFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const fvMesh& mesh = patch().boundaryMesh().mesh();
    const fvPatch& patch = this->patch();

    if (mesh.foundObject<volTensorField>(kappaName_))
    {
        const fvPatchField<tensor>& k = 
            patch.lookupPatchField<volTensorField, tensor>
        (
            kappaName_
        );
        
        vectorField n = (patch.Sf() / patch.magSf());
        
        gradient() = -q_ / ((n & k) & n);
    }
    else
    {
        const fvPatchField<scalar>& k = 
            patch.lookupPatchField<volScalarField, scalar>
        (
            kappaName_
        );

        gradient() = -q_ / k;
    }

    fvPatchField<scalar>::updateCoeffs();
}


void fixedHeatFluxFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);

#ifdef OPENFOAMFOUNDATION
    writeEntryIfDifferent<word>(os, "kappa", "k", kappaName_);
    writeEntry(os, "heatFlux", q_);
    writeEntry(os, "gradient", gradient());
    writeEntry(os, "value", *this);
#elif OPENFOAMESI
    os.writeEntryIfDifferent<word>("kappa", "k", kappaName_);
    q_.writeEntry("heatFlux", os);
    gradient().writeEntry("gradient", os);
    this->writeEntry("value", os);
#endif   
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, fixedHeatFluxFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
