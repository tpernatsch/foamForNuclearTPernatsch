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
    Foam::plenumSpringPressureFvPatchVectorField

Description
    The `plenumSpringPressure` boundary condition replicates the counteracting 
    force of the plenum spring.

    Physically, the total pressure applied on the patch is calculated as the sum 
    of the gap pressure (derived from a `gapGasModel`) and the reaction force of the 
    spring. The spring force is determined based on the spring modulus and the 
    relative displacement of the surfaces to which the spring is attached. 

    The shear stress on the patch is set to zero.

    \note
    This BC is designed to be applied to the top of the pellet column (i.e., the 
    topmost fuel patch) and, if present, the inner surface of the top cladding cap. 
    For models without a top cap (where the cladding ends with an open annular 
    cylindrical surface), the user can use the `topCladRingPressure` BC to account 
    for the spring effect.

    \warning
    This BC assumes the presence of a `gapGasModel`.

________________________________________________________________________________

\par Options
    The `plenumSpringPressure` fvPatchField can be selected in the patch subdictionary 
    inside the `boundaryField` subdictionary of the displacement field.

    Parameters in the patch subdictionary for `plenumSpringPressure`:

    - <b>`springModulus`</b> - Modulus of the plenum spring in N/m, representing the 
    stiffness of the spring.
    - <b>`initialSpringLoading`</b> - Pre-loading of the plenum spring in meters.
    - <b>`fuelTopPatches`</b> - List of patches attached to the bottom end of the spring 
    (typically the top surface of the fuel column).
    - <b>`topCapInnerPatches`</b> - List of patches attached to the top end of the spring 
    (either the inner surface of the top cap, if present, or the top annular surface 
    of the cladding in an open cladding model).

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
    To apply the `plenumSpringPressure` boundary condition to a given patch, you 
    need to specify its typename in the patch definition in the boundary condition 
    file (e.g., `0/D` or `0/DD` for simulations involving incremental mechanical 
    solvers).

    Example: Applying `plenumSpringPressure` to the `fuelTop` patch:

    \verbatim
        fuelTop
        {
            type            plenumSpringPressure;

            springModulus           1e4;
            initialSpringLoading    0;

            fuelTopPatches          (fuelTop);
            topCapInnerPatches      (cladTop);

            // Activate plane strain approximation (default is false)
            planeStrain     false;

            // Relaxation factor
            relax           1.0;

            // Initial displacement value
            value           $internalField;
        }
    \endverbatim

SourceFiles
    plenumSpringPressureFvPatchVectorField.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
    Switzerland, Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef plenumSpringPressureFvPatchVectorField_H
#define plenumSpringPressureFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "tractionDisplacementFvPatchVectorField.H"
#include "plenumSpringBase.H"
#include "gapGasModel.H" 


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class plenumSpringPressureFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class plenumSpringPressureFvPatchVectorField
:
    public tractionDisplacementFvPatchVectorField,
    public plenumSpringBase
{

    // Private Data

public:   

    //- Runtime type information
    TypeName("plenumSpringPressure");


    // Constructors

        //- Construct from patch and internal field
        plenumSpringPressureFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        plenumSpringPressureFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  plenumSpringPressureFvPatchVectorField onto a new patch
        plenumSpringPressureFvPatchVectorField
        (
            const plenumSpringPressureFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );


    // Member functions

        // Access Functions

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
