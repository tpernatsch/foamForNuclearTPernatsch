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
    The `coolantPressure` boundary condition is used to apply the coolant pressure 
    (fixed or time-dependent) to a patch. It is a simplified re-implementation of 
    the `tractionDisplacement` boundary condition, allowing the user to work 
    directly with the coolant pressure instead of specifying both traction and 
    pressure values.

    \note
    This BC inherits from the parten class `tractionDisplacement` the possibility
    of activating a plane strain approximation as well as a spring-dashpot system.

________________________________________________________________________________

\par Options
    The `coolantPressure` fvPatchField can be selected in the patch subdictionary 
    inside the `boundaryField` subdictionary of the displacement field.

    Parameters in the patch subdictionary for `coolantPressure`:

    - <b>`coolantPressureList`</b> - A time-dependent pressure table. Each entry 
    consists of a time step (using the time unit defined in `solverDict`) followed 
    by the pressure value in Pa.
    - <b>`coolantPressure`</b> - Fixed pressure in Pa. Used only if 
    `coolantPressureList` is not provided.
    - <b>`planeStrain`</b> - Activates the plane strain approximation for the normal 
    stress at the boundary. When enabled, the normal strain is assumed constant 
    across the last layer of cells. **Default: `false`.**
    - <b>`fixedSpring`</b> - Activates a fixed spring-dashpot system for additional 
    stability. **Default: `false`.**
    - <b>`fixedSpringModulus`</b> - Spring modulus in N/m. Required when `fixedSpring` 
    is set to `true`.
    - <b>`dashpotModulus`</b> - Dashpot modulus in N/m. Required when `fixedSpring` 
    is set to `true`.
    - <b>`relax`</b> - Relaxation factor for gradient updates. **Default: `1.0`.**
    - <b>`value`</b> - Initial displacement value (not stress).

________________________________________________________________________________

Usage
    To apply the `coolantPressure` boundary condition to a given patch, you need 
    to specify its typename in the patch definition in the boundary condition file 
    (e.g., `0/D` or `0/DD` for simulations involving incremental mechanical 
    solvers).

    Example: Applying `coolantPressure` to the `cladOuter` patch:

    \verbatim
        cladOuter
        {
            type            coolantPressure;

            // Fixed pressure in case the coolantPressureList is not provided
            // coolantPressure 1e5;

            coolantPressureList
            {
                type                table;

                // Use "file" or directly provide "values"
                file                "$FOAM_CASE/constant/lists/coolantPressureList";
                
                // Optional: values can be provided directly
                // values
                // (
                //     (0.0 1e5)
                //     (100.0 1e6)
                // );

                format              foam;      // Optional data format
                outOfBounds         clamp;     // Optional handling for out-of-bounds values
                interpolationScheme linear;    // Optional interpolation method
            }

            // Activate plane strain approximation (default is false)
            planeStrain     false;

            // Relaxation factor
            relax           1.0;

            // Initial displacement value
            value           $internalField;
        }
    \endverbatim

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

        //- Update the pressure and traction
        virtual void updateTraction();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
