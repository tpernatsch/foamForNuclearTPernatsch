/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "codedTractionDisplacementFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

codedTractionDisplacementFvPatchVectorField::
codedTractionDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    codedFixedGradientFvPatchVectorField(p, iF),
    traction_(p.size(), vector::zero),
    relax_(1),
    stressName_("sigma"),
    pressureList_()
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = vector::zero;
}


codedTractionDisplacementFvPatchVectorField::
codedTractionDisplacementFvPatchVectorField
(
    const codedTractionDisplacementFvPatchVectorField& tdpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    codedFixedGradientFvPatchVectorField(tdpvf, p, iF, mapper),
    traction_(tdpvf.traction_),
    relax_(tdpvf.relax_),
    stressName_(tdpvf.stressName_),
    pressureList_()
{}


codedTractionDisplacementFvPatchVectorField::
codedTractionDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    codedFixedGradientFvPatchVectorField(p, iF, dict),        
    traction_(p.size(), vector::zero),
    relax_(dict.lookupOrDefault<scalar>("relax", 1)),
    stressName_(dict.lookupOrDefault<word>("stressName", "sigma")),    
    pressureList_()
{ 
    if (dict.found("value"))
    {
        fvPatchVectorField::operator=
        (
            vectorField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchVectorField::operator=(patchInternalField());
    }
    // Info << "if " <<     patch().lookupPatchField<volVectorField, vector>(internalField().name()) << endl;
    gradient() = vector::zero;
    // this->evaluate();
}


codedTractionDisplacementFvPatchVectorField::
codedTractionDisplacementFvPatchVectorField
(
    const codedTractionDisplacementFvPatchVectorField& tdpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    codedFixedGradientFvPatchVectorField(tdpvf, iF),
    traction_(tdpvf.traction_),
    relax_(tdpvf.relax_),
    stressName_(tdpvf.stressName_),
    pressureList_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void codedTractionDisplacementFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchVectorField::autoMap(m);
    m(traction_, traction_);
}


void codedTractionDisplacementFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchVectorField::rmap(ptf, addr);

    const codedTractionDisplacementFvPatchVectorField& dmptf =
        refCast<const codedTractionDisplacementFvPatchVectorField>(ptf);

    traction_.rmap(dmptf.traction_, addr);
}

void codedTractionDisplacementFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Make sure library containing user-defined fvPatchField is up-to-date
    updateLibrary(name_);

    const fvPatchVectorField& fvp = redirectPatchField();

    const_cast<fvPatchVectorField&>(fvp).updateCoeffs();
    
    // Copy the traction from the on-the fly bc through value
    traction_=(fvp);

    const fvPatchField<scalar>& mu =
        patch().lookupPatchField<volScalarField, scalar>("mu");

    const fvPatchField<scalar>& lambda =
        patch().lookupPatchField<volScalarField, scalar>("lambda");
        
    const fvPatchField<symmTensor>& stress =
        patch().lookupPatchField<volSymmTensorField, symmTensor>(stressName_);

    // Lookup the gradient field
    const fvPatchField<tensor>& gradField =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad" + internalField().name()
        );        
    
    scalarField twoMuLambda(2*mu + lambda);
    
    vectorField n(patch().nf());  

    gradient() = 
    (
        (
              relax_*(traction_ - (n & stress)) 
        ) / twoMuLambda + (n &gradField)
    );
    
    fixedGradientFvPatchVectorField::updateCoeffs();
}

void codedTractionDisplacementFvPatchVectorField::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    // Lookup the gradient field
    const fvPatchField<tensor>& gradField =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad" + internalField().name()
        );

    // Face unit normal
    const vectorField n = patch().nf();

    // NOTE:: Because of how OF calculates delta, the correction implemented 
    // here is useless (vector k is by construct equal to zero).
    
    // Delta vectors
    const vectorField delta = patch().delta();

    // Non-orthogonal correction vectors
    vectorField k = delta - n*(n&delta);

    Field<vector>::operator=
        (
            patchInternalField()
            + (k & gradField.patchInternalField())
            + gradient()/patch().deltaCoeffs()
        );

    fvPatchField<vector>::evaluate();
}


void codedTractionDisplacementFvPatchVectorField::write(Ostream& os) const
{
    codedFixedGradientFvPatchVectorField::write(os);
    writeEntry(os, "value". traction_);
    writeEntryIfDifferent<scalar>(os, "relax", 1, relax_);
    writeEntryIfDifferent<word>(os, "stressName", "sigma", stressName_);
    writeEntry(os, "value");
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    codedTractionDisplacementFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
