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

#include "radiativeConvectiveSinkFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

radiativeConvectiveSinkFvPatchScalarField::radiativeConvectiveSinkFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvPatchField<scalar>(p, iF),
    H_(p.size(), 0.0),
    emissivity_(p.size(), 0.0),
    T0_(p.size(), 0.0),
    kappaName_("k"),
    alpha_(p.size(), 0.0)
{}


radiativeConvectiveSinkFvPatchScalarField::radiativeConvectiveSinkFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<scalar>(p, iF),
    H_("alpha", dict, p.size()),
    emissivity_("emissivity", dict, p.size()),
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


radiativeConvectiveSinkFvPatchScalarField::radiativeConvectiveSinkFvPatchScalarField
(
    const radiativeConvectiveSinkFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper,
    const bool mappingRequired
)
:
    fvPatchField<scalar>(ptf, p, iF, mapper),
    H_(p.size()),
    emissivity_(p.size()),
    T0_(p.size()),
    kappaName_(ptf.kappaName_),
    alpha_(p.size())
{
     if (mappingRequired)
     {
         if (mapper.hasUnmapped())
         {
             H_ = 0.0;
             emissivity_ = 0.0;
             T0_ = 0.0;
             alpha_ = 0.0;
         }
         
    #ifdef OPENFOAMFOUNDATION    
        mapper(H_, H_);
        mapper(emissivity_, emissivity_);
        mapper(T0_, T0_);
        mapper(alpha_, alpha_);
    #elif OPENFOAMESI
        H_.autoMap(mapper);
        emissivity_.autoMap(mapper);
        T0_.autoMap(mapper);
        alpha_.autoMap(mapper);
    #endif
    }
}


radiativeConvectiveSinkFvPatchScalarField::radiativeConvectiveSinkFvPatchScalarField
(
    const radiativeConvectiveSinkFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvPatchField<scalar>(tppsf, iF),
    H_(tppsf.H_),
    emissivity_(tppsf.emissivity_),
    T0_(tppsf.T0_),
    kappaName_(tppsf.kappaName_),
    alpha_(tppsf.alpha_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void radiativeConvectiveSinkFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fvPatchField<scalar>::autoMap(m);
#ifdef OPENFOAMFOUNDATION    
    m(H_, H_);
    m(emissivity_, emissivity_);
    m(T0_, T0_);
    m(alpha_, alpha_);
#elif OPENFOAMESI
    H_.autoMap(m);
    emissivity_.autoMap(m);
    T0_.autoMap(m);
    alpha_.autoMap(m);
#endif    

}


void radiativeConvectiveSinkFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fvPatchField<scalar>::rmap(ptf, addr);

    const radiativeConvectiveSinkFvPatchScalarField& dmptf =
        refCast<const radiativeConvectiveSinkFvPatchScalarField>(ptf);

    H_.rmap(dmptf.H_, addr);
    emissivity_.rmap(dmptf.emissivity_, addr);
    T0_.rmap(dmptf.T0_, addr);
    alpha_.rmap(dmptf.alpha_, addr);
}


void radiativeConvectiveSinkFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const fvMesh& mesh = patch().boundaryMesh().mesh();
    const scalarField& Tb = *this;
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
        scalarField alphaN = (H_ + 5.67e-8*emissivity_*(Tb*Tb + T0_*T0_)*(Tb+T0_));
        
        alpha_ = alphaP / (alphaP + alphaN);
    }
    else
    {
        const fvPatchField<scalar>& k = patch.lookupPatchField<volScalarField, scalar>
        (
            kappaName_
        );

        scalarField alphaP = k*deltaCoeffs;
        scalarField alphaN = (H_ + 5.67e-8*emissivity_*(Tb*Tb + T0_*T0_)*(Tb+T0_));
        
        alpha_ = alphaP / (alphaP + alphaN);
    }

    fvPatchField<scalar>::updateCoeffs();
}


void radiativeConvectiveSinkFvPatchScalarField::evaluate
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
radiativeConvectiveSinkFvPatchScalarField::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<scalarField>(new scalarField(alpha_));
}


tmp<Field<scalar> >
radiativeConvectiveSinkFvPatchScalarField::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return (1-alpha_)*T0_;
}


tmp<Field<scalar> >
radiativeConvectiveSinkFvPatchScalarField::gradientInternalCoeffs() const
{
    const scalarField& deltaCoeffs = patch().deltaCoeffs();

    return (alpha_-1)*deltaCoeffs;
}


tmp<Field<scalar> >
radiativeConvectiveSinkFvPatchScalarField::gradientBoundaryCoeffs() const
{
    const scalarField& deltaCoeffs = patch().deltaCoeffs();

    return (1-alpha_)*deltaCoeffs*T0_;
}


void radiativeConvectiveSinkFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
#ifdef OPENFOAMFOUNDATION
    writeEntryIfDifferent<word>(os, "kappa", "k", kappaName_);
    writeEntry(os, "T0", T0_);
    writeEntry(os, "alpha", H_);
    writeEntry(os, "emissivity", emissivity_);
    writeEntry(os, "value", *this);
#elif OPENFOAMESI
    os.writeEntryIfDifferent<word>("kappa", "k", kappaName_);
    T0_.writeEntry("T0", os);
    H_.writeEntry("alpha", os);
    emissivity_.writeEntry("emissivity", os);
    this->writeEntry("value", os);
#endif
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, radiativeConvectiveSinkFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
