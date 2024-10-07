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
    Foam::timeDependentAxialProfileHTCfvPatchScalarField

Description
    Heat sink boundary defined by an axial profile of convective heat transfer 
    coefficients and a sink temperature axial profile. 

    \htmlinclude fvPatchField_timeDependentAxialProfileHTCfvPatchScalarField.html

SourceFiles
    timeDependentAxialProfileHTCfvPatchScalarField.C

\mainauthor
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\contribution
    A. Scolaro, E. Brunetto, C. Fiorina - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
    Switzerland, Laboratory for Reactor Physics and Systems Behaviour)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef timeDependentAxialProfileHTCfvPatchScalarField_H
#define timeDependentAxialProfileHTCfvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "axialProfile.H"
#include "InterpolateTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class timeDependentAxialProfileHTCfvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class timeDependentAxialProfileHTCfvPatchScalarField
:
    public fvPatchField<scalar>
{
    // Private data

        //- Tabulated axial locations
        scalarField zValues_;
        
        //- Tabulated time values
        const scalarField timeValues_;
        
        //- Tabulated HTC Profile data
        FieldField<Field, scalar> hProfileData_;

        //- Tabulated Coolant Profile data
        FieldField<Field, scalar> T0ProfileData_;
    
        //- Axial location interpolation method
        interpolateTableBase::interpolationMethod zMethod_;
        
        //- Time interpolation method
        interpolateTableBase::interpolationMethod tMethod_;
        
        //- h Interpolation table
        scalarFieldInterpolateTable hTable_;

        //- T0 Interpolation table
        scalarFieldInterpolateTable T0Table_;

        //- Axial profile dict
        dictionary dict_;

        //- Name of thermal conductivity field
        word kappaName_;
        
        //- T0 patch profile @ current timestep - working variable
        scalarField T0_;

        //- h patch profile @ current timestep - for write()
        scalarField h_;
        
        //- Boundary heat transfer factor - working variable
        scalarField alpha_;
        

public:

    //- Runtime type information
    TypeName("timeDependentAxialProfileHTC");


    // Constructors

        //- Construct from patch and internal field
        timeDependentAxialProfileHTCfvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        timeDependentAxialProfileHTCfvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&,
            const bool valueRequired=true
        );

        //- Construct by mapping given timeDependentAxialProfileHTCfvPatchScalarField
        //  onto a new patch
        timeDependentAxialProfileHTCfvPatchScalarField
        (
            const timeDependentAxialProfileHTCfvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&,
            const bool mappingRequired=true
        );

        //- Construct as copy setting internal field reference
        timeDependentAxialProfileHTCfvPatchScalarField
        (
            const timeDependentAxialProfileHTCfvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new timeDependentAxialProfileHTCfvPatchScalarField(*this, iF)
            );
        }


    // Member functions
    
        //- Return true if this patch field fixes a value.
        //  Needed to check if a level has to be specified while solving
        //  Poissons equations.
        virtual bool fixesValue() const
        {
            return true;
        }

        // Mapping functions
            
        
        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Evaluate the patch field
        virtual void evaluate
        (
            const Pstream::commsTypes commsType=Pstream::commsTypes::blocking
        );

        //- Return the matrix diagonal coefficients corresponding to the evaluation of the value of this patchField
        virtual tmp<Field<scalar> > valueInternalCoeffs
        (
            const tmp<scalarField>&
        ) const;

        //- Return the matrix source coefficients corresponding to the evaluation of the value of this patchField,
        virtual tmp<Field<scalar> > valueBoundaryCoeffs
        (
            const tmp<scalarField>&
        ) const;

        //- Return the matrix diagonal coefficients corresponding to the evaluation of the gradient of this patchField
        virtual tmp<Field<scalar> > gradientInternalCoeffs() const;

        //- Return the matrix source coefficients corresponding to the evaluation of the gradient of this patchField
        virtual tmp<Field<scalar> > gradientBoundaryCoeffs() const;
        
        //- Write to an output stream
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
