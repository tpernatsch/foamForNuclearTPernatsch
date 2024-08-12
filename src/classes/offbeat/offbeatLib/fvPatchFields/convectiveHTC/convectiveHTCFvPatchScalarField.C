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

#include "convectiveHTCFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

convectiveHTCFvPatchScalarField::convectiveHTCFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvPatchField<scalar>(p, iF),
    h_(p.size(), 0.0),
    T0_(p.size(), 0.0),
    kappaName_("k"),
    alpha_(p.size(), 0.0)
{}


convectiveHTCFvPatchScalarField::convectiveHTCFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<scalar>(p, iF),
    h_("h", dict, p.size()),
    T0_("T0", dict, p.size()),
    kappaName_(dict.lookupOrDefault<word>("kappa", "k")),
    alpha_(p.size(), 0.0)
{
    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
}


convectiveHTCFvPatchScalarField::convectiveHTCFvPatchScalarField
(
    const convectiveHTCFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchField<scalar>(ptf, p, iF, mapper),
    h_(ptf.h_),
    T0_(ptf.T0_),
    kappaName_(ptf.kappaName_),
    alpha_(ptf.alpha_)
{}


convectiveHTCFvPatchScalarField::convectiveHTCFvPatchScalarField
(
    const convectiveHTCFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvPatchField<scalar>(tppsf, iF),
    h_(tppsf.h_),
    T0_(tppsf.T0_),
    kappaName_(tppsf.kappaName_),
    alpha_(tppsf.alpha_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void convectiveHTCFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fvPatchField<scalar>::autoMap(m);
#ifdef OPENFOAMFOUNDATION
    m(h_, h_);
    m(T0_, T0_);
    m(alpha_, alpha_);
#elif OPENFOAMESI    
    h_.autoMap(m);
    T0_.autoMap(m);
    alpha_.autoMap(m);
#endif 
}


void convectiveHTCFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fvPatchField<scalar>::rmap(ptf, addr);

    const convectiveHTCFvPatchScalarField& dmptf =
        refCast<const convectiveHTCFvPatchScalarField>(ptf);

    h_.rmap(dmptf.h_, addr);
    T0_.rmap(dmptf.T0_, addr);
    alpha_.rmap(dmptf.alpha_, addr);
}


void convectiveHTCFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const fvMesh& mesh = patch().boundaryMesh().mesh();
    const fvPatch& patch = this->patch();
    const scalarField& deltaCoeffs = patch.deltaCoeffs();

    if (mesh.foundObject<volTensorField>(kappaName_))
    {
        const fvPatchField<tensor>& k = patch.lookupPatchField<volTensorField, tensor>
        (
            kappaName_
        );
        
        vectorField n = (patch.Sf() / patch.magSf());
        scalarField alphaP = ((n & k) & n)*deltaCoeffs;
        const scalarField& alphaN = h_;
        
        alpha_ = alphaP / (alphaP + alphaN);
    }
    else
    {
        const fvPatchField<scalar>& k = patch.lookupPatchField<volScalarField, scalar>
        (
            kappaName_
        );

        scalarField alphaP = k*deltaCoeffs;
        const scalarField& alphaN = h_;
        
        alpha_ = alphaP / (alphaP + alphaN);
    }

    fvPatchField<scalar>::updateCoeffs();
}


void convectiveHTCFvPatchScalarField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }
    
    fvPatchField<scalar>::operator=(
            alpha_*patchInternalField() + (1-alpha_)*T0_
    );
    
    fvPatchField<scalar>::evaluate();
}


tmp<Field<scalar> >
convectiveHTCFvPatchScalarField::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<scalarField>(new scalarField(alpha_));
}


tmp<Field<scalar> >
convectiveHTCFvPatchScalarField::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return (1-alpha_)*T0_;
}


tmp<Field<scalar> >
convectiveHTCFvPatchScalarField::gradientInternalCoeffs() const
{
    const scalarField& deltaCoeffs = patch().deltaCoeffs();

    return (alpha_-1)*deltaCoeffs;
}


tmp<Field<scalar> >
convectiveHTCFvPatchScalarField::gradientBoundaryCoeffs() const
{
    const scalarField& deltaCoeffs = patch().deltaCoeffs();

    return (1-alpha_)*deltaCoeffs*T0_;
}


void convectiveHTCFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
#ifdef OPENFOAMFOUNDATION
    writeEntryIfDifferent<word>(os, "kappa", "k", kappaName_);
    writeEntry(os, "h", h_);
    writeEntry(os, "T0", T0_);
    writeEntry(os, "value", *this);
#elif OPENFOAMESI
    os.writeEntryIfDifferent<word>("kappa", "k", kappaName_);
    os.writeEntry("h", h_);
    os.writeEntry("T0", T0_);
    os.writeEntry("value", *this);
#endif    
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, convectiveHTCFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
