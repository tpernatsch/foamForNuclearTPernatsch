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
    Foam::gapContactFvPatchVectorField

Description
    Daugther-class to the tractionDisplacement BC for modeling the contact between 
    fuel and cladding. It uses the penalty method and a Coulomb friction model.

    \htmlinclude fvPatchField_gapContactFvPatchVectorField.html

SourceFiles
    gapContactFvPatchVectorField.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
    Switzerland, Laboratory for Reactor Physics and Systems Behaviour)\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\contribution
    E. Brunetto, C. Fiorina - EPFL

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef gapContactFvPatchVectorField_H
#define gapContactFvPatchVectorField_H

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

class gapContactFvPatchVectorField
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
            inline const gapContactFvPatchVectorField& neighbour() const
            {
                return refCast<const gapContactFvPatchVectorField>
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
    TypeName("gapContact");


    // Constructors

        //- Construct from patch and internal field
        gapContactFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        gapContactFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  gapContactFvPatchVectorField onto a new patch
        gapContactFvPatchVectorField
        (
            const gapContactFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        gapContactFvPatchVectorField
        (
            const gapContactFvPatchVectorField&,
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
                new gapContactFvPatchVectorField(*this, iF)
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
