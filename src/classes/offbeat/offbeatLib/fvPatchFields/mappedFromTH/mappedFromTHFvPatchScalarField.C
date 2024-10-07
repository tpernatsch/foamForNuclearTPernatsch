 /*---------------------------------------------------------------------------*\
   =========                 |
   \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
    \\    /   O peration     |
     \\  /    A nd           | www.openfoam.com
      \\/     M anipulation  |
 -------------------------------------------------------------------------------
     Copyright (C) 2020-2021 OpenCFD Ltd.
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
 
 #include "mappedFromTHFvPatchScalarField.H"
 #include "volFields.H"
 #include "interpolationCell.H"
 #include "addToRunTimeSelectionTable.H"

 
 // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
 

 Foam::mappedFromTHFvPatchScalarField::mappedFromTHFvPatchScalarField
 (
     const fvPatch& p,
     const DimensionedField<scalar, volMesh>& iF
 )
 :
     mixedFvPatchField<scalar>(p, iF)
 {
     this->refValue() = Zero;
     this->refGrad() = Zero;
     this->valueFraction() = 0.0;
 }
 
 

 Foam::mappedFromTHFvPatchScalarField::mappedFromTHFvPatchScalarField
 (
     const fvPatch& p,
     const DimensionedField<scalar, volMesh>& iF,
     const dictionary& dict
 )
 :
     // Bypass dictionary constructor (all reading handled later)
     // but cannot use NO_READ since will still trigger an evaluate()
     mixedFvPatchField<scalar>(p, iF)
 {
     fvPatchFieldBase::readDict(dict);  // Consistent with a dict constructor
 
     this->readValueEntry(dict, IOobjectOption::MUST_READ);
 
     if (this->readMixedEntries(dict))
     {
         // Full restart
     }
     else
     {
         // Start from user entered data. Assume fixedValue.
         this->refValue() = *this;
         this->refGrad() = Zero;
         this->valueFraction() = 1.0;
     }
 
 // This blocks (crashes) with more than two worlds!
 //
 }
 
 
 Foam::mappedFromTHFvPatchScalarField::mappedFromTHFvPatchScalarField
 (
     const mappedFromTHFvPatchScalarField& ptf,
     const fvPatch& p,
     const DimensionedField<scalar, volMesh>& iF,
     const fvPatchFieldMapper& mapper
 )
 :
     mixedFvPatchField<scalar>(ptf, p, iF, mapper)
 
{}
 Foam::mappedFromTHFvPatchScalarField::mappedFromTHFvPatchScalarField
 (
     const mappedFromTHFvPatchScalarField& ptf
 )
 :
     mixedFvPatchField<scalar>(ptf)
 {}
 

 Foam::mappedFromTHFvPatchScalarField::mappedFromTHFvPatchScalarField
 (
     const mappedFromTHFvPatchScalarField& ptf,
     const DimensionedField<scalar, volMesh>& iF
 )
 :
     mixedFvPatchField<scalar>(ptf, iF)
 {}
 
 
 // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
 
 void Foam::mappedFromTHFvPatchScalarField::autoMap(const fvPatchFieldMapper& m)
 {
     mixedFvPatchField<scalar>::autoMap(m);
 }
 

 void Foam::mappedFromTHFvPatchScalarField::rmap
 (
     const fvPatchField<scalar>& ptf,
     const labelList& addr
 )
 {
     mixedFvPatchField<scalar>::rmap(ptf, addr);
 }


 void Foam::mappedFromTHFvPatchScalarField::updateCoeffs()
 {
     if (this->updated())
     {
         return;
     }
 
     // Set my weights

     scalarField tk(patch().lookupPatchField<volScalarField, Foam::scalar>("k"));

     scalarField myWeights(tk*patch().deltaCoeffs());


     //Set nbrWeight

     scalarField nbrWeights(patch().lookupPatchField<volScalarField, Foam::scalar>("mappedhtc"));
 
     this->refValue() = patch().lookupPatchField<volScalarField, Foam::scalar>("mappedTFluid");
     this->refGrad() = Zero;
     this->valueFraction() = nbrWeights/(nbrWeights + myWeights);
 
     mixedFvPatchField<scalar>::updateCoeffs();
 
 }
 
 
 void Foam::mappedFromTHFvPatchScalarField::write(Ostream& os) const
 {
     mixedFvPatchField<scalar>::write(os);
 }
 
 
 // ************************************************************************* //

  // * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(mappedFromTHFvPatchScalarField, 0);
addToPatchFieldRunTimeSelection
(
    fvPatchScalarField,
    mappedFromTHFvPatchScalarField
);
}

// ************************************************************************* //