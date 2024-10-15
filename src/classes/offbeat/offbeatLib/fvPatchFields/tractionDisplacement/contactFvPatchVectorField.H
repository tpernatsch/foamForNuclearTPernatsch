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
    Applies interface pressure if surface penetration is detected. Requires regioncoupled.
    If there is no penetration, it applies a constant pressure. 

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
        
        //- Penalty factor
        scalar penaltyScaleFact_;
        
        //- Gap width [m]
        mutable scalarField gapWidth_;

        //- Interface pressure
        mutable scalarField interfaceP_;

        //- Outer pressure
        mutable scalarField outerPressure_;

        //- Apply force on slave only if concentric?
        Switch forceConcentricSlave;

        //- Relax InterfacePressure
        scalar relaxInterP_;

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
