/*---------------------------------------------------------------------------*\
   =========                 |
   \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
    \\    /   O peration     | Website:  https://openfoam.org
     \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
 
 Class
     Foam::porosityOutletFvPatchScalarField
 
 Description
     Boundary condition implemented for testing porosity outflow. From 
     'advective' BC.
 
 SourceFiles
     porosityOutletFvPatchScalarField.C
 
 \*---------------------------------------------------------------------------*/
 
 #ifndef porosityOutletFvPatchScalarField_H
 #define porosityOutletFvPatchScalarField_H
 
 #include "mixedFvPatchFields.H"
 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
 namespace Foam
 {
 
 /*---------------------------------------------------------------------------*\
                     Class porosityOutletFvPatchScalarField Declaration
 \*---------------------------------------------------------------------------*/
 
 class porosityOutletFvPatchScalarField
 :
    public mixedFvPatchScalarField
 {
 protected:
 
    // Private Data

    scalar factor_;
 
 
 public:
 
    //- Runtime scalar information
    TypeName("porosityOutlet");

 
     // Constructors
 
        //- Construct from patch and internal field
        porosityOutletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        porosityOutletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given porosityOutletFvPatchScalarField
        //  onto a new patch
        porosityOutletFvPatchScalarField
        (
            const porosityOutletFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        porosityOutletFvPatchScalarField(const porosityOutletFvPatchScalarField&) = delete;

        //- Copy constructor setting internal field reference
        porosityOutletFvPatchScalarField
        (
            const porosityOutletFvPatchScalarField&,
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
                new porosityOutletFvPatchScalarField(*this, iF)
            );
        }
 
 
     // Member Functions
 
         // Access
 
 
         // Evaluation functions
 
            //- Calculate and return the advection speed at the boundary
            virtual scalarField advectionSpeed() const;

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

 
        //- Write
        virtual void write(Ostream&) const;
 };
 
 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
 } // End namespace Foam
 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 

 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
 #endif
 
 // ************************************************************************* //
