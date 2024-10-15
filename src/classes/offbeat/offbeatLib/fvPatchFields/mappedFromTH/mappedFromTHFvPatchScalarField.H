 /*---------------------------------------------------------------------------*\
   =========                 |
   \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
    \\    /   O peration     |
     \\  /    A nd           | www.openfoam.com
      \\/     M anipulation  |
 -------------------------------------------------------------------------------
     Copyright (C) 2019-2020 OpenCFD Ltd.
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
 
 Class
     Foam::mappedFromTHFvPatchScalarField
 
 Group
     grpGenericBoundaryConditions grpCoupledBoundaryConditions
 
 Description
     This boundary condition maps the value at a set of cells or patch faces
     back to *this.
 
     The sample mode is set by the underlying mapping engine, provided by the
     mappedPatchBase class.
 
 Usage
     \table
         Property     | Description                      | Required | Default
         field        | Name of field to be mapped       | no  | this field name
         weightField  | Name of field to use as weight   | no |
         setAverage   | Use average value                | no  | false |
         average      | Average value to use if \c setAverage = yes | partly |
     \endtable
 
     Example of the boundary condition specification:
     \verbatim
     <patchName>
     {
         type            mappedMixed;
         field           T;              // optional field name
 
         value           uniform 273.0;
 
         refValue        $value;
         refGradient     uniform 0.0;
         valueFraction   uniform 1.0;
     }
     \endverbatim
 
 Note
     Supports multi-world operation:
 
     +----+                      +----+
     |    |                      |    |
     |    |                      |    |
     +----+                      +----+
     worldA                      worldB
     regionA                     regionB
     patchA                      patchB
 
 
 See also
     Foam::mappedPatchBase
     Foam::mappedPolyPatch
     Foam::mappedFvPatch
     Foam::mappedFieldFvPatchField
     Foam::mixedFvPatchField
 
 SourceFiles
     mappedFromTHFvPatchScalarField.C
 
 \*---------------------------------------------------------------------------*/
 
 #ifndef mappedFromTHFvPatchScalarField_H
 #define mappedFromTHFvPatchScalarField_H
 
 #include "mixedFvPatchFields.H"
 #include "scalarField.H"
 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
 namespace Foam
 {
 
 /*---------------------------------------------------------------------------*\
                   Class mappedFromTHFvPatchScalarField Declaration
 \*---------------------------------------------------------------------------*/
 
 class mappedFromTHFvPatchScalarField
 :
     public mixedFvPatchField<scalar>
 {
 

 
 public:
 
     //- Runtime type information
     TypeName("mappedFromTH");
 
 
     // Constructors
 
         //- Construct from patch and internal field
         mappedFromTHFvPatchScalarField
         (
             const fvPatch&,
             const DimensionedField<scalar, volMesh>&
         );
 
         //- Construct from patch, internal field and dictionary
         mappedFromTHFvPatchScalarField
         (
             const fvPatch&,
             const DimensionedField<scalar, volMesh>&,
             const dictionary&
         );
 
         //- Construct by mapping given
         //  mappedFromTHFvPatchScalarField
         //  onto a new patch
         mappedFromTHFvPatchScalarField
         (
             const mappedFromTHFvPatchScalarField&,
             const fvPatch&,
             const DimensionedField<scalar, volMesh>&,
             const fvPatchFieldMapper&
         );
 
         //- Construct as copy
         mappedFromTHFvPatchScalarField
         (
             const mappedFromTHFvPatchScalarField&
         );
 
         //- Construct and return a clone
         virtual tmp<fvPatchField<scalar>> clone() const
         {
             return tmp<fvPatchField<scalar>>
             (
                 new mappedFromTHFvPatchScalarField
                 (
                     *this
                 )
             );
         }
 
         //- Construct as copy setting internal field reference
         mappedFromTHFvPatchScalarField
         (
             const mappedFromTHFvPatchScalarField&,
             const DimensionedField<scalar, volMesh>&
         );
 
         //- Construct and return a clone setting internal field reference
         virtual tmp<fvPatchField<scalar>> clone
         (
             const DimensionedField<scalar, volMesh>& iF
         ) const
         {
             return tmp<fvPatchField<scalar>>
             (
                 new mappedFromTHFvPatchScalarField
                 (
                     *this,
                     iF
                 )
             );
         }
 
 
     // Member functions
 
         // Mapping functions
 
             //- Map (and resize as needed) from self given a mapping object
             virtual void autoMap
             (
                 const fvPatchFieldMapper&
             );
 
             //- Reverse map the given fvPatchField onto this fvPatchField
             virtual void rmap
             (
                 const fvPatchField<scalar>&,
                 const labelList&
             );


 
 
         // Evaluation functions
 
             //- Update the coefficients associated with the patch field
             virtual void updateCoeffs();
 
 
         //- Write
         virtual void write(Ostream&) const;
 };
 
 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
 } // End namespace Foam
 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
 /*#ifdef NoRepository
     #include "mappedFromTHFvPatchScalarField.C"
 #endif*/
 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
 #endif
 
 // ************************************************************************* //