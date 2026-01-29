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
    Foam::gapPressureFvPatchVectorField

Description
    The `gapPressure` boundary condition applies the gap pressure to the selected 
    patch. It is derived from the `tractionDisplacement` boundary condition and 
    shares its general mechanism (including the spring-dashpot system). 
    However, with `gapPressure`, the user does not need to specify the traction 
    vector or pressure directly.

    \warning
    This BC is specifically designed to work in conjunction with a `gapGasModel`
    that dynamically provides the gap pressure at every iteration.

________________________________________________________________________________

\par Options
    The `gapPressure` fvPatchField can be selected in the patch subdictionary 
    inside the `boundaryField` subdictionary of the displacement field.

    Parameters in the patch subdictionary for `gapPressure`:

    - <b>`planeStrain`</b> - Activates the plane strain approximation for the 
    normal stress at the boundary. When enabled, the normal strain is assumed 
    constant across the last layer of cells. **Default: `false`.**
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
    To apply the `gapPressure` boundary condition to a given patch, you need to 
    specify its typename in the patch definition in the boundary condition file 
    (e.g., `0/D` or `0/DD` for simulations involving incremental mechanical 
    solvers).

    Example: Applying `gapPressure` to the `fuelInner` patch:

    \verbatim
        fuelInner
        {
            type            gapPressure;

            // Activate plane strain approximation for the normal stress at the
            // boundary. It is false by default.
            planeStrain     false;

            // Relaxation factor for the gradient update
            relax           1.0;

            // Initial displacement value
            value           $internalField;
        }
    \endverbatim

SourceFiles
    gapPressureFvPatchVectorField.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
    Switzerland, Laboratory for Reactor Physics and Systems Behaviour)\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\contribution
    E. Brunetto, C. Fiorina - EPFL

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef gapPressureFvPatchVectorField_H
#define gapPressureFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "tractionDisplacementFvPatchVectorField.H"
#include "gapGasModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class gapPressureFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class gapPressureFvPatchVectorField
:
    public tractionDisplacementFvPatchVectorField
{

    // Private Data
  

public:   

    //- Runtime type information
    TypeName("gapPressure");


    // Constructors

        //- Construct from patch and internal field
        gapPressureFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        gapPressureFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  gapPressureFvPatchVectorField onto a new patch
        gapPressureFvPatchVectorField
        (
            const gapPressureFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        gapPressureFvPatchVectorField
        (
            const gapPressureFvPatchVectorField&,
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
                new gapPressureFvPatchVectorField(*this, iF)
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