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
    Foam::unilateralContactFvPatchVectorField

Description
    The `unilateralContact` boundary condition enforces a contact-like behavior 
    on the patch. In the parallel direction, sliding is always allowed, with no 
    shear forces applied. 

    In the normal direction, the patch displacement is restricted from being 
    positive (relative to the patch normal). If a positive displacement is attempted, 
    the displacement is forced to be zero in the normal direction.
    The user can specify a relaxation factor that can be used to aid 
    convergence. The spring constraint is enforced in an implicit manner.

    When the patch detaches inward (negative displacement relative to the patch normal), 
    the traction/pressure is applied. This traction/pressure can be read 
    from the patchField dictionary. Otherwise it is assumed to be zero

________________________________________________________________________________

\par Options
    The `unilateralContact` fvPatchField can be selected in the patch subdictionary 
    inside the `boundaryField` subdictionary of the displacement field.

    Parameters in the patch subdictionary for `unilateralContact`:
    
    - <b>`traction`</b> - Fixed traction vector in Pa. Zero by default.
    - <b>`pressure`</b> - Fixed pressure in Pa. Zero by default.
    - <b>`relax`</b> - Relaxation factor for the contact contribution. Helps improve 
      convergence. **Default: `1`.**

________________________________________________________________________________

Usage
    To apply the `unilateralContact` boundary condition to a given patch, the 
    following example can be used as a template:

    \verbatim
        cladOuter
        {
            type            unilateralContact;

            // Relaxation factor for the contact contribution
            relax           1;

            // Pressure for the part of the boundary not in contact (zero by default)
            pressure uniform 0.0;

            // Traction for the part of the boundary not in contact (zero by default)
            traction uniform (0 0 0);
        }
    \endverbatim

SourceFiles
    unilateralContactFvPatchVectorField.C

\mainauthor
    Alessandro Scolaro, Milica Krstovic - EPFL


\*---------------------------------------------------------------------------*/

#ifndef unilateralContactFvPatchVectorField_H
#define unilateralContactFvPatchVectorField_H

#include "fvPatchFields.H"
#include "solidDirectionMixedFvPatchVectorField.H"
#include "Switch.H"
#include "Table.H"
#include "gapGasModel.H"

#ifdef OPENFOAMESI
    #include "Function1.H"
#endif   

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class unilateralContactFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class unilateralContactFvPatchVectorField
:
    public solidDirectionMixedFvPatchVectorField
{
protected:
    
    // Protected data
  
        //- Surface traction vector
        vectorField traction_;
        
        //- Surface pressure
        scalarField pressure_;
  
        //- Total traction vector (pressure + traction)
        vectorField totalTraction_;
        
        //- Old total traction vector (pressure + traction)
        //  Only for incremental solvers
        autoPtr<vectorField> totalTraction0_;
        
        // relaxation factor for update
        scalar relax_;

        //- Spring modulus
        scalar fixedSpringModulus_; 

        // //- Patch weights for spring-dashpot system
        // scalarField w_;   

        //- Value of current time step (in s)
        scalar currentTime_;

        //- Provide the normal to the boundary plane (the direction along which
        // there is no allowed penetration through)
        vector planeNormal_;

        //- Provide the point through which the plane passes
        vector planePoint_;

        //- Relative (to cell size) width for linear ramp of contact stiffness
        scalar contactRelRampWidth_;

        //- Penetration relative (to cell size) tolerance for contact plane
        scalar contactRelOffset_;

        //- Penalty scale factor for contact stiffness
        scalar penaltyScaleFactor_;

public:

    //- Runtime type information
    TypeName("unilateralContact");


    // Constructors

        //- Construct from patch and internal field
        unilateralContactFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        unilateralContactFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&,
            const bool valueRequired = true
        );

        //- Construct by mapping given
        // unilateralContactFvPatchVectorField
        //  onto a new patch
        unilateralContactFvPatchVectorField
        (
            const unilateralContactFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&,
            const bool mappingRequired=true
        );

        //- Construct as copy setting internal field reference
        unilateralContactFvPatchVectorField
        (
            const unilateralContactFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<vector> > clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<vector> >
            (
                new unilateralContactFvPatchVectorField(*this, iF)
            );
        }

    // Destructor

        virtual ~unilateralContactFvPatchVectorField()
        {}


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
                const fvPatchField<vector>&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;


    // Member operators

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
