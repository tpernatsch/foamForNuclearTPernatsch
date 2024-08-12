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
    Foam::fixedTorqueFvPatchVectorField

Description
    Fixed torque in Nm uniformly distributed over the boundary faces.

SourceFiles
    fixedTorqueFvPatchVectorField.C

\mainauthor
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    January 2022

\*---------------------------------------------------------------------------*/

#ifndef fixedTorqueFvPatchVectorField_H
#define fixedTorqueFvPatchVectorField_H

#include "tractionDisplacementFvPatchVectorField.H"
#include "quaternion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class fixedTorqueFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class fixedTorqueFvPatchVectorField
:
    public tractionDisplacementFvPatchVectorField
{

    // Private Data
    
        //- User-specified torque
        scalar torque_;
        
        //- Axis of rotation
        vector axis_;
        
        //- Handle rotation around axis
        quaternion rotation_;
        
        //- Centre of rotation
        vector C0_;
  
    // Private Member Functions
        
        //- Check that the provided axis is valid
        void check();
        
        //- Calculation the centre of rotation (centroid of boundary faces)
        void calculateRotationCentre();

public:   

    //- Runtime type information
    TypeName("fixedTorque");


    // Constructors

        //- Construct from patch and internal field
        fixedTorqueFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedTorqueFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  fixedTorqueFvPatchVectorField onto a new patch
        fixedTorqueFvPatchVectorField
        (
            const fixedTorqueFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        fixedTorqueFvPatchVectorField
        (
            const fixedTorqueFvPatchVectorField&,
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
                new fixedTorqueFvPatchVectorField(*this, iF)
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
                const fvPatchVectorField&,
                const labelList&
            );


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
