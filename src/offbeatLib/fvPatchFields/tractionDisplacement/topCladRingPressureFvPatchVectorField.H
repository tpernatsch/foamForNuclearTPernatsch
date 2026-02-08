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
    Foam::topCladRingPressureFvPatchVectorField

Description
    The `topCladRingPressure` boundary condition is used for the top cladding 
    annular ring when the top cap is not explicitly modeled (often done even in 
    2D or 3D full rod simulations, either due to the difficulty of finding 
    precise cap dimensions or to save computational time).

    Physically, this boundary condition ensures equilibrium by balancing the forces 
    acting on the annular ring at the top of the cladding tube. If one imagines to
    remove the top cap, the forces acting on the remaining cladding ring (internal and 
    external pressures, as well as the spring force) must balance for the tube to remain 
    in equilibrium. The `topCladRingPressure` BC approximates this force distribution 
    as uniform across the ring surface, using the following contributions:

    - **Outer fluid pressure** (downward force): Specified as a fixed value or a 
    time-dependent list.
    - **Gap pressure** (upward force): Acting on the plenum side of the gap.
    - **Plenum spring pressure** (upward force).

    \note
    In the real rod, the coolant and the gap gas and plenum spring act on 
    different surfaces: the coolant force acts on the outer top surface of the cap, 
    while the gap gas and spring forces act on the inner bottom surface of the cap. 
    To ensure accuracy, these pressure contributions are scaled by the ratio of the 
    corresponding surface areas before being applied uniformly to the top cladding 
    ring surface. For this reason, the user must provide the dimensions of the top cap 
    (internal and external radii).

    \note
    This BC inherits from the parten class `tractionDisplacement` the possibility
    of activating a plane strain approximation as well as a spring-dashpot system.

    \warning
    This BC must be used in combination to a `plenumSpringPressure` BC on another patch,
    typically the fuel top most patch. In the respective dictionary for this patch,
    `topCapInnerPatches` list must reference the patch using `topCladRingPressure`.

________________________________________________________________________________

\par Options
    The `topCladRingPressure` fvPatchField can be selected in the patch subdictionary 
    inside the `boundaryField` subdictionary of the displacement field.

    Parameters in the patch subdictionary for `topCladRingPressure`:

    - <b>`topCapRadius`</b> - Radius of the top cap on the fluid side (outer radius) in meters.
    - <b>`topCapPlenumRadius`</b> - Radius of the top cap on the plenum side (inner radius) in meters.
    - <b>`coolantPressureList`</b> - A time-dependent coolant pressure table. Each entry 
    consists of a time step (using the time unit defined in `solverDict`) followed 
    by the pressure value in Pa.
    - <b>`coolantPressure`</b> - Fixed coolant pressure in Pa. Used only if `coolantPressureList` 
    is not provided.- <b>`tractionList`</b> - A time-dependent traction table. Each entry consists of 
    a time step and a corresponding traction vector value in Pa.

    Parameters in the patch subdictionary derived from `tractionDisplacement`:
    
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
    To apply the `topCladRingPressure` boundary condition to a given patch, you need 
    to specify its typename in the patch definition in the boundary condition file 
    (e.g., `0/D` or `0/DD` for simulations involving incremental mechanical solvers).

    Example 1: Applying `topCladRingPressure` to the `cladTop` patch

    \verbatim
        cladTop
        {
            type            topCladRingPressure;

            topCapRadius       0.00525;
            topCapPlenumRadius 0.0045;

            // Fixed coolant pressure if coolantPressureList is not provided
            // coolantPressure 1e5;

            coolantPressureList
            {
                type                table;

                // Use "file" or directly provide "values"
                file                "$FOAM_CASE/constant/lists/topCladdingRingPressureCoolantList";
                
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

        // Accompanying `plenumSpringPressure` BC definition for `fuelTop` patch
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
    topCladRingPressureFvPatchVectorField.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
    Switzerland, Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef topCladRingPressureFvPatchVectorField_H
#define topCladRingPressureFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "coolantPressureFvPatchVectorField.H"
#include "plenumSpringBase.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class topCladRingPressureFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class topCladRingPressureFvPatchVectorField
:
    public coolantPressureFvPatchVectorField
{

    // Private Data

        //- Top cap radius on the fluid side
        scalar topCapRadius_;

        //- Top cap radius on the plenum side
        scalar topCapPlenumRadius_;

        //- Label of fuel top patch
        label fuelTopPatchID_;

        //- Name of strain field
        word displacementName_;   

public:   

    //- Runtime type information
    TypeName("topCladRingPressure");


    // Constructors

        //- Construct from patch and internal field
        topCladRingPressureFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        topCladRingPressureFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  topCladRingPressureFvPatchVectorField onto a new patch
        topCladRingPressureFvPatchVectorField
        (
            const topCladRingPressureFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );


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
