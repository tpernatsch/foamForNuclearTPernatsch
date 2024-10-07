/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::tractionDisplacementFvPatchVectorField

Description
    Fixed traction boundary condition for the standard linear elastic, fixed
    coefficient displacement equation.

    \htmlinclude fvPatchField_tractionDisplacementFvPatchVectorField.html

SourceFiles
    tractionDisplacementFvPatchVectorField.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
    Switzerland, Laboratory for Reactor Physics and Systems Behaviour)\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\contribution
    E. Brunetto, C. Fiorina - EPFL

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef tractionDisplacementFvPatchVectorField_H
#define tractionDisplacementFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "Table.H"

#ifdef OPENFOAMESI
    #include "Function1.H"
#endif    


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class tractionDisplacementFvPatch Declaration
\*---------------------------------------------------------------------------*/

class tractionDisplacementFvPatchVectorField
:
    public fixedGradientFvPatchVectorField
{

    // Private Data
        
protected:

    // Protected Data
  
        // Surface traction vector
        vectorField traction_;
        
        // Surface pressure
        scalarField pressure_;
        
        // relaxation factor for update
        scalar relax_;

        //- Name of stress field
        word stressName_;   

        //- Apply plane strain correction
        Switch planeStrain_;

        //- Value of current time step (in s)
        scalar currentTime_;

#ifdef OPENFOAMFOUNDATION
        //- Time varying list of pressure values
        autoPtr<Function1s::Table<scalar> > pressureList_;

        //- Time varying list of traction values
        autoPtr<Function1s::Table<vector> > tractionList_;

#elif OPENFOAMESI
        typedef Foam::Function1Types::Table<scalar> scalarTable;
        typedef Foam::Function1Types::Table<vector> vectorTable;

        //- Time varying list of pressure values
        autoPtr<scalarTable > pressureList_;

        //- Time varying list of traction values
        autoPtr<vectorTable >  tractionList_;
#endif        

        //- Apply fixed spring-dashpot system or not
        Switch fixedSpring_;

        //- Spring modulus
        scalar fixedSpringModulus_; 

        //- Dashpot modulus
        scalar dashpotModulus_; 

        //- Patch weights for spring-dashpot system
        scalarField w_;
        
public:

    //- Runtime type information
    TypeName("tractionDisplacement");


    // Constructors

        //- Construct from patch and internal field
        tractionDisplacementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary.
        //  The bool is used to decide whether to initialize the fields
        tractionDisplacementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&,
            const bool valueRequired = true
        );

        //- Construct by mapping given
        //  tractionDisplacementFvPatchVectorField onto a new patch
        tractionDisplacementFvPatchVectorField
        (
            const tractionDisplacementFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        tractionDisplacementFvPatchVectorField
        (
            const tractionDisplacementFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new tractionDisplacementFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            virtual const vectorField& traction() const
            {
                return traction_;
            }

            virtual vectorField& traction()
            {
                return traction_;
            }

            virtual const scalarField& pressure() const
            {
                return pressure_;
            }

            virtual  scalarField& pressure()
            {
                return pressure_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchVectorField&,
                const labelList&
            );


        // Evaluation functions

            //- Return gradient at boundary
            virtual tmp<Field<vector>> snGrad() const;

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Evaluate the patch
            virtual void evaluate(const Pstream::commsTypes);
 
            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<vector>> gradientInternalCoeffs() const;
 
            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<vector>> gradientBoundaryCoeffs() const; 
 
            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<vector>> valueInternalCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<vector>> valueBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
