/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::coolantPressureFvPatchVectorField

Description
    Re-implementation of the tractionDisplacement mother class that allows the user
    to deal only with "coolantPressure" instead than the more confusion "traction" and
    "pressure". 

    \htmlinclude fvPatchField_coolantPressureFvPatchVectorField.html 

SourceFiles
    coolantPressureFvPatchVectorField.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
    Switzerland, Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef coolantPressureFvPatchVectorField_H
#define coolantPressureFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "tractionDisplacementFvPatchVectorField.H"

#ifdef OPENFOAMESI
    #include "Function1.H"
#endif     

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class coolantPressureFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class coolantPressureFvPatchVectorField
:
    public tractionDisplacementFvPatchVectorField
{

protected:

    // Protected Data

        //- Fluid pressure
        scalarField coolantPressure_;

#ifdef OPENFOAMFOUNDATION
        //- Fluid pressure time-dependent list
        autoPtr<Function1s::Table<scalar>> coolantPressureList_;
#elif OPENFOAMESI
        typedef Foam::Function1Types::Table<scalar> scalarTable;
        
        //- Fluid pressure time-dependent list
        autoPtr<scalarTable> coolantPressureList_;
#endif          

public:   

    //- Runtime type information
    TypeName("coolantPressure");


    // Constructors

        //- Construct from patch and internal field
        coolantPressureFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        coolantPressureFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  coolantPressureFvPatchVectorField onto a new patch
        coolantPressureFvPatchVectorField
        (
            const coolantPressureFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        coolantPressureFvPatchVectorField
        (
            const coolantPressureFvPatchVectorField&,
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
                new coolantPressureFvPatchVectorField(*this, iF)
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
