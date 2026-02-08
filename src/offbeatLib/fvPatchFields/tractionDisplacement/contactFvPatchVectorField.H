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
    Foam::contactFvPatchVectorField

Description
    The `contact` boundary condition uses the penalty method to handle
    contact between two coupled surfaces of type `regionCoupledOFFBEAT`. This 
    boundary condition is designed to manage both normal and tangential (frictional) 
    forces between contacting surfaces. 

    \note
    This fvPatchField is originally partly based on the [solids4Foam](https://www.solids4foam.com/) 
    library.
    
    \warning
    This boundary condition assumes that the same boundary condition is applied to a
    pair of `regionCoupledOFFBEAT` patches, which are  coupled using the AMI
    (Arbitrary Mesh Interface) mapping algorithm. It is strongly suggested
    to assign the `owner` role to the smallest surface, e.g. the `fuelOuter` surface
    in most cases.
    
    For the normal component, if no penetration is detected, the pressure on both sides 
    of the gap equals the gap pressure. If penetration occurs, a contact pressure 
    proportional to the penetration is added, calculated as:
        
    $$
    \begin{aligned}
    ip = -f \cdot K \cdot g
    \end{aligned}
    $$

    where:

    - \f$ip\f$: Interface normal pressure.
    - \f$g\f$: Penetration or `gapWidth` (negative if penetration occurs).
    - \f$K\f$: the minimum surface stiffness between those of master and slave faces.
    On each side, the stiffness is calculated internally as \f$K = K_m \cdot A/V\f$,
    where \f$K_m\f$ is the bulk modulus, and \f$A\f$ and \f$V\f$ are the face area and 
    volume of the cell.
    - \f$f\f$: Penalty scale factor, introduced to facilitate convergence at 
    the cost of allowing larger penetration.

    \note
    The penalty scale factor \f$f\f$ should be chosen as large as possible (typically 
    close to 1) while ensuring convergence.

    For the tangential (friction) component a Coulomb model is used. The 
    magnitude of the friction vector is limited by the normal interface pressure
    times a friction coefficient and it is given by:
        
    $$
    \begin{aligned}
    |\tau| = \min(\alpha \cdot ip, |f_{fr} \cdot \mu \cdot \textbf{s}|)
    \end{aligned}
    $$

    where:

    - \f$\tau\f$: The frictional vector component.
    - \f$\alpha\f$: Friction coefficient.
    - \f$\mu\f$: Minimum shear stiffness between those of slave and master face.
    On each side the shear stiffness is calculated internally as \f$\mu = \mu_m \cdot A/V\f$, 
    where \f$\mu_m\f$ is the shear modulus, and \f$A\f$ and \f$V\f$ are the face area 
    and cell volume.
    - \f$s\f$: Slip vector, the cumulative tangential displacement between the two 
    surfaces.
    - \f$f_{fr}\f$: Penalty scale factor for friction, introduced to improve 
    convergence at the cost of larger penetration.

    The friction vector is oriented opposite to the slip direction and added 
    to the normal pressure component to compute the total surface traction vector.

    \note
    The gap width calculation considers the displacement of face centers (`D` for 
    total mechanics solvers, `DD` for incremental solvers) so it correctly calculates
    the updated gap size even if a `smallStrain` mechanics solver is used
    
    \note 
    This BC allows the user to introduce an `offset` to adjust the gap size for 
    sensitivity analysis.

________________________________________________________________________________

\par Options
    The `contact` fvPatchField can be selected in the patch subdictionary
    inside the `boundaryField` subdictionary of the displacement field.

    Parameters in the patch subdictionary for `contact`:

    - <b>`patchType`</b> - Specifies the patch field type. **Must be set to 
    `regionCoupledOFFBEAT`.**
    - <b>`penaltyFactor`</b> - Penalty scale factor for the normal pressure. 
    Controls the convergence and penetration balance. **Default: `0.1`.**
    Typical values: `0.01` to `1`.
    - <b>`penaltyFactorFriction`</b> - Penalty scale factor for the frictional 
    forces. **Default: `0.1`.** Typical values: `0.01` to `1`.
    - <b>`frictionCoefficient`</b> - Coefficient of friction between surfaces. 
    **Default: `0`.**
    - <b>`relaxInterfacePressure`</b> - Relaxation factor for the normal pressure 
    update. **Default: `1`.** Typical values: `0.01` to `1`.
    - <b>`relaxFriction`</b> - Relaxation factor for the frictional force update. 
    **Default: `1`.** Typical values: `0.01` to `1`.
    - <b>`offset`</b> - Offset added to the gap width calculation. A positive value 
    increases the gap size. **Default: `0`.**
    - <b>`rigidMasterNormal`</b> - If `true`, the normal component is applied only 
    on the slave side in rigid master mode. **Default: `false`.**
    - <b>`rigidMasterFriction`</b> - If `true`, the friction component is applied 
    only on the slave side in rigid master mode. **Default: `false`.**

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
    To apply the `contact` boundary condition to a given patch, specify its
    typename in the patch definition in the boundary condition file (e.g., `0/D` or 
    `0/DD` for incremental mechanical solvers).

    Example: Applying `contact` to the `fuelOuter` and `cladInner` patches:

    \verbatim
        "fuelOuter|cladInner"
        {
            type            contact;
            patchType       regionCoupledOFFBEAT;

            penaltyFactor             0.01;
            penaltyFactorFriction     0.01;

            frictionCoefficient        0.1;

            relaxInterfacePressure     0.01;
            relaxFriction              0.01;

            // Offset for gap width calculation (a positive value increases the gap size)
            offset          0.0;

            // If rigidMasterNormal(Friction) is true, the normal (friction) component
            // is applied only on the slave side
            rigidMasterNormal       false;
            rigidMasterFriction     false;

            // Activate plane strain approximation for the normal stress at the boundary
            planeStrain     false;

            // Relaxation factor for the gradient update
            relax           1.0;

            // Initial displacement value
            value           $internalField;
        }
    \endverbatim

SourceFiles
    contactFvPatchVectorField.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
    Switzerland, Laboratory for Reactor Physics and Systems Behaviour)\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\contribution
    E. Brunetto, C. Fiorina - EPFL

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef contactFvPatchVectorField_H
#define contactFvPatchVectorField_H

#include "tractionDisplacementFvPatchVectorField.H"
#include "regionCoupledBaseOFFBEATFvPatch.H"
#include "gapGasModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class tractionDisplacementFvPatch Declaration
\*---------------------------------------------------------------------------*/

class contactFvPatchVectorField
:
    public tractionDisplacementFvPatchVectorField
{

    // Private Data

        //- Local reference to region couple patch
        const regionCoupledBaseOFFBEATFvPatch& regionCoupledPatch_;
        
        //- Penalty scale factor
        scalar penaltyScaleFact_;
        
        //- Penalty scale factor for friction
        scalar penaltyScaleFactFric_;   

        //- Penetration tolerance for augmentedLagrangian approach
        scalar penetrationTolerance_;      
        
        //- Coulomb friction coefficient
        scalar frictionCoeff_;
        
        //- Gap width [m]
        mutable scalarField gapWidth_;

        //- Interface pressure
        mutable scalarField interfaceP_;

        //- Slip between contact patches
        mutable vectorField slip_;

        //- Old value of slip between contact patches
        mutable vectorField slip0_;

        //- Friction component of contact traction
        mutable vectorField frictionTraction_;

        //- Relax InterfacePressure
        scalar relaxInterP_;

        //- Relax friction force
        scalar relaxFriction_;

        //- Takes track of current time value
        scalar currentTime_;

        //- Offset for gapWidth (a positive value increases the gap width)
        scalar offset_;

        //- If true, normal traction applies only on slave patch (typically cladding)
        Switch rigidMasterNormal_;

        //- If true, tangetial traction applies only on slave patch (typically cladding)
        Switch rigidMasterFriction_;


    // Private functions
        
protected:

    // Protected functions
    
            //- Return the neighbour patchField
            inline const contactFvPatchVectorField& neighbour() const
            {
                return refCast<const contactFvPatchVectorField>
                (
                    regionCoupledPatch_.neighbFvPatch().lookupPatchField
                    <volVectorField, vector>
                    (
                        internalField().name()
                    )
                );
            }

            //- Return the gap width
            tmp<scalarField> gapWidth() const;

            //- Minimum patch stiffness for interface pressure evaluation
            scalarField boundaryStiffness() const;

            //- Minimum patch shear stiffness for friction evaluation
            scalarField boundaryShearStiffness() const;
            
public:

    //- Runtime type information
    TypeName("contact");


    // Constructors

        //- Construct from patch and internal field
        contactFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        contactFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  contactFvPatchVectorField onto a new patch
        contactFvPatchVectorField
        (
            const contactFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        contactFvPatchVectorField
        (
            const contactFvPatchVectorField&,
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
                new contactFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access
        
            //- Interface pressure
            inline const scalarField& interfaceP() const
            {
                return interfaceP_;
            }      
        
            //- Interface pressure
            inline const vectorField& frictionTraction() const
            {
                return frictionTraction_;
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
            
        // Modify

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
