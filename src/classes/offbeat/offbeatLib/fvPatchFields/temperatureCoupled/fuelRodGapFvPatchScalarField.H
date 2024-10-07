/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
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
    Foam::fuelRodGapFvPatchScalarField

Description
    Region coupled patchField for temperature fields in solid heat conduction.

Usage
    \verbatim
    fuelOuter
    {
        type            resistiveGap;
        patchType       regionCoupledOFFBEAT;
        coupled         true;
        roughness       uniform 1e-6;
        relax           1.0;
        relaxGapWidth   1.0;
        value           $internalField;
    }
    \endverbatim

SourceFiles
    fuelRodGapFvPatchScalarField.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
    Switzerland, Laboratory for Reactor Physics and Systems Behaviour)\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\contribution
    E. Brunetto, C. Fiorina - EPFL

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef fuelRodGapFvPatchScalarField_H
#define fuelRodGapFvPatchScalarField_H

#include "resistiveGapFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class fuelRodGapFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class fuelRodGapFvPatchScalarField
:
    public resistiveGapFvPatchScalarField
{
    // Private data
    
        //- Gap width [m]
        mutable scalarField gapWidth_;

        //- Gap conductance
        mutable scalarField hGap_;
        
        //- Surface roughness [m]
        scalarField R_;
        
        //- Name of emissivity field
        word emissivityName_;

        //- Name of gap gas model
        word gapGasModelName_;    

        //- Relaxation factor for hGap
        scalar relax_;    

        //- Relaxation factor for gapWidth
        scalar relaxGapWidth_;    

        
    // Private functions

protected:

    // Protected functions
    
            //- Return the neighbour patchField
            inline const fuelRodGapFvPatchScalarField& neighbour() const
            {
                return refCast<const fuelRodGapFvPatchScalarField>
                (
                    regionCoupledPatch_.neighbFvPatch().lookupPatchField
                    <volScalarField, scalar>
                    (
                        internalField().name()
                    )
                );
            }

            //- Return the gap width
            tmp<scalarField> gapWidth() const;
            
            //- Return the interface conductance
            tmp<scalarField> hGap() const;

public:

    //- Runtime type information
    TypeName("fuelRodGap");


    // Constructors

        //- Construct from patch and internal field
        fuelRodGapFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fuelRodGapFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fuelRodGapFvPatchScalarField
        // onto a new patch
        fuelRodGapFvPatchScalarField
        (
            const fuelRodGapFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        fuelRodGapFvPatchScalarField
        (
            const fuelRodGapFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<scalar> > clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<scalar> >
            (
                new fuelRodGapFvPatchScalarField(*this, iF)
            );
        }


    //- Destructor
    virtual ~fuelRodGapFvPatchScalarField()
    {}


    // Member functions

        // Access
        
            //- Return the surface roughness
            inline const scalarField& R() const
            {
                return R_;
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
                const fvPatchScalarField&,
                const labelList&
            );
            
        // Evaluation functions

            //- Update the coefficients associated with the patch field
            //  Sets Updated to true
            virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
