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

#include "zeroStressGradientFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

zeroStressGradientFvPatchVectorField::
zeroStressGradientFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    tractionDisplacementFvPatchVectorField(p, iF)
{}


zeroStressGradientFvPatchVectorField::
zeroStressGradientFvPatchVectorField
(
    const zeroStressGradientFvPatchVectorField& tdpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    tractionDisplacementFvPatchVectorField(tdpvf, p, iF, mapper)
{}


zeroStressGradientFvPatchVectorField::
zeroStressGradientFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict,
    const bool valueRequired
)
:
    tractionDisplacementFvPatchVectorField(p, iF, dict, false)     
{}    


zeroStressGradientFvPatchVectorField::
zeroStressGradientFvPatchVectorField
(
    const zeroStressGradientFvPatchVectorField& tdpvf,
    const DimensionedField<vector, volMesh>& iF
)
:  
    tractionDisplacementFvPatchVectorField(tdpvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void zeroStressGradientFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchVectorField::autoMap(m);
}


void zeroStressGradientFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchVectorField::rmap(ptf, addr);
}


void zeroStressGradientFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    vectorField n(patch().nf());
    vectorField nCurrent(n);
    scalarField magSf(patch().magSf());

    if(patch().boundaryMesh().mesh().foundObject<volTensorField>("Finv"))
    {
        const tensorField& Finv =
        patch().boundaryMesh().mesh().foundObject<fvMesh>("referenceMesh") ?
        patch().lookupPatchField<volTensorField, tensor>("relFinv") :
        patch().lookupPatchField<volTensorField, tensor>("Finv");
        
        // Patch unit normals (deformed configuration)
        nCurrent = (Finv.T() & n);
        nCurrent /= mag(nCurrent);
    };

    const fvPatchField<scalar>& mu =
        patch().lookupPatchField<volScalarField, scalar>("mu");

    const fvPatchField<scalar>& lambda =
        patch().lookupPatchField<volScalarField, scalar>("lambda");

    // Lookup the explicit component of the stress
    const fvPatchField<tensor>& sigmaExp =
        patch().lookupPatchField<volTensorField, tensor>("sigmaExp");

    symmTensorField stressI
    ( 
        db().lookupObject<volSymmTensorField>(stressName_).internalField(), 
        patch().faceCells()
    );
    
    scalarField twoMuLambda(2*mu + lambda);

    vectorField newGradient
    (
        (
            (n & stressI) - (n & sigmaExp)
        ) / twoMuLambda
    );

    // Relax gradient
    gradient() = (1-relax_)*gradient() + relax_*newGradient;

    fixedGradientFvPatchVectorField::updateCoeffs();
}


void zeroStressGradientFvPatchVectorField::write(Ostream& os) const
{
    fixedGradientFvPatchVectorField::write(os);

#ifdef OPENFOAMFOUNDATION
    writeEntryIfDifferent<scalar>(os, "relax", 1, relax_);
    writeEntryIfDifferent<word>(os, "stressName", "sigma", stressName_);
    writeEntry(os, "value", *this);
#elif OPENFOAMESI
    os.writeEntryIfDifferent<scalar>("relax", 1, relax_);
    os.writeEntryIfDifferent<word>( "stressName", "sigma", stressName_);
    this->writeEntry("value", os);
#endif    
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    zeroStressGradientFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
