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
    The `fixedTorque` boundary condition applies a fixed torque (in Nm) uniformly 
    distributed over the boundary faces. This boundary condition is useful for 
    simulating rotational forces acting on a boundary.

    The applied torque is defined by the user in the boundary dictionary, specifying 
    the magnitude, axis of rotation, and optionally, the center of rotation. If no 
    center is provided, the rotation center is automatically calculated.

________________________________________________________________________________

\par Options
    The `fixedTorque` fvPatchField can be selected in the patch subdictionary 
    inside the `boundaryField` subdictionary of the displacement field.

    Parameters in the patch subdictionary for `fixedTorque`:

    - <b>`torque`</b> - Magnitude of the applied torque in Nm.
    - <b>`axis`</b> - Vector specifying the axis of rotation (e.g., `(0 0 1)` for 
    the z-axis).
    - <b>`centre`</b> - (Optional) Vector specifying the center of rotation in meters. 
    If not provided, the rotation center is automatically calculated.    

    Parameters in the patch subdictionary derived from `tractionDisplacement`:

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
    To apply the `fixedTorque` boundary condition to a given patch, you need to 
    specify its typename in the patch definition in the boundary condition file 
    (e.g., `0/D` or `0/DD` for simulations involving incremental mechanical 
    solvers).

    Example: Applying `fixedTorque` to the `rotatingBoundary` patch:

    \verbatim
        rotatingBoundary
        {
            type            fixedTorque;

            // Magnitude of the applied torque in Nm
            torque          10.0;

            // Axis of rotation
            axis            (0 0 1);

            // (Optional) Center of rotation
            centre          (0 0 0);

            // Activate plane strain approximation (default is false)
            planeStrain     false;

            // Relaxation factor for the gradient update
            relax           1.0;

            // Initial displacement value
            value           $internalField;
        }
    \endverbatim

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
