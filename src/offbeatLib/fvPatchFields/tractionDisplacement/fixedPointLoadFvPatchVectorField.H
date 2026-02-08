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
    Foam::fixedPointLoadFvPatchVectorField

Description
    The `fixedPointLoad` boundary condition applies a fixed point load to a patch. 
    The specified force is distributed over the boundary faces closest to the 
    user-defined point. This BC is useful for simulating localized forces acting 
    on specific regions of a boundary, trying to replicate FEM boundary conditions
    that apply the force on a single node.

    The user specifies the location of the point and the magnitude of the force 
    in the boundary dictionary. If the force magnitude is set to zero, the 
    simulation will terminate with an error.

________________________________________________________________________________

\par Options
    The `fixedPointLoad` fvPatchField can be selected in the patch subdictionary 
    inside the `boundaryField` subdictionary of the displacement field.

    Parameters in the patch subdictionary for `fixedPointLoad`:

    - <b>`point`</b> - Vector specifying the location of the point in meters.
    - <b>`force`</b> - Vector specifying the magnitude and direction of the applied 
    force in N. The force is distributed over the faces closest to the specified 
    point.
    - <b>`pressure`</b> - (Optional) Initial pressure field for restart simulations.

    Parameters in the patch subdictionary derived from `tractionDisplacement`:
	
    - <b>`fixedSpring`</b> - Activates a fixed spring-dashpot system for additional 
    stability. **Default: `false`.**
    - <b>`fixedSpringModulus`</b> - Spring modulus in N/m. Required when `fixedSpring` 
    is set to `true`.
    - <b>`dashpotModulus`</b> - Dashpot modulus in N/m. Required when `fixedSpring` 
    is set to `true`.
	- <b>`planeStrain`</b> - Activates the plane strain approximation for the normal 
	stress at the boundary. When enabled, the normal strain is assumed constant 
	across the last layer of cells. **Default: `false`.**
	- <b>`relax`</b> - Relaxation factor for gradient updates. **Default: `1.0`.**
	- <b>`value`</b> - Initial displacement value (not stress).

________________________________________________________________________________

Usage
    To apply the `fixedPointLoad` boundary condition to a given patch, you need 
    to specify its typename in the patch definition in the boundary condition 
    file (e.g., `0/D` or `0/DD` for simulations involving incremental mechanical 
    solvers).

    Example: Applying `fixedPointLoad` to the `loadBoundary` patch:

    \verbatim
        loadBoundary
        {
            type            fixedPointLoad;

            // Location of the point in meters
            point           (0.5 0.5 0);

            // Magnitude and direction of the applied force in N
            force           (0 0 -100);

            // (Optional) Initial pressure field for restart simulations
            pressure        uniform 0;

            // Activate plane strain approximation (default is false)
            planeStrain     false;

            // Relaxation factor for the gradient update
            relax           1.0;

            // Initial displacement value
            value           $internalField;
        }
    \endverbatim

SourceFiles
    fixedPointLoadFvPatchVectorField.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
    Switzerland, Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef fixedPointLoadFvPatchVectorField_H
#define fixedPointLoadFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "tractionDisplacementFvPatchVectorField.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class fixedPointLoadFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class fixedPointLoadFvPatchVectorField
:
    public tractionDisplacementFvPatchVectorField
{

    // Private Data

        //- Point where the load is applied
        vector point_;

        //- Force
        vector force_;    

public:   

    //- Runtime type information
    TypeName("fixedPointLoad");


    // Constructors

        //- Construct from patch and internal field
        fixedPointLoadFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedPointLoadFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  fixedPointLoadFvPatchVectorField onto a new patch
        fixedPointLoadFvPatchVectorField
        (
            const fixedPointLoadFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        fixedPointLoadFvPatchVectorField
        (
            const fixedPointLoadFvPatchVectorField&,
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
                new fixedPointLoadFvPatchVectorField(*this, iF)
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
