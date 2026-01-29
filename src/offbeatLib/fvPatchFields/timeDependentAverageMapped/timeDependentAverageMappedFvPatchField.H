/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::timeDependentAverageMappedFvPatchField

Description
    mappedFixedValue boundary conditions where the average is defined using
    a time-dependent interpolation table.

Usage
    \table
        Property     | Description             | Required    | Default value
        fieldName    | name of field to be mapped | no       | this field name
        timePoints   | time points for interpolation table (scalarField) | yes | 
        averagePoints | average value correspoding to each time point (Field<Type>) | yes |
        interpolationScheme | type of interpolation scheme | no |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            timeDependentAverageMapped;
        fieldName       T;
        interpolationScheme cell;
        timePoints      ( 0 1 2 3 ... );
        averagePoints   ( 1 1 1.5 ... );
        value           uniform 0;
    }
    \endverbatim

See also
    Foam::timeDependentAverageMappedFvPatchField

SourceFiles
    timeDependentAverageMappedFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef timeDependentAverageMappedFvPatchField_H
#define timeDependentAverageMappedFvPatchField_H

#include "fixedValueFvPatchFields.H"
#include "mappedPatchFieldBase.H"
#include "InterpolateTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class timeDependentAverageMappedFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class timeDependentAverageMappedFvPatchField
:
    public fixedValueFvPatchField<Type>,
    public mappedPatchFieldBase<Type>
{

protected:

    // Protected data
    
        typedef InterpolateTable<Field<Type>, Type, Type> 
            interpolateTable;
        
        typedef interpolateTableBase::interpolationMethod interpolationMethod;
    
        //- List of time points
        scalarField timePoints_;
        
        //- List of average values
        Field<Type> averageValues_;

        //- Method for interpolation
        interpolationMethod method_;

        //- Time-dependent table for rod-average LHGR (W/m)
        interpolateTable table_;

    // Protected Member Functions

        const mappedPatchBase& mapper
        (
            const fvPatch& p,
            const DimensionedField<Type, volMesh>& iF
        );

public:

    //- Runtime type information
    TypeName("timeDependentMapped");

    // Constructors

        //- Construct from patch and internal field
        timeDependentAverageMappedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        timeDependentAverageMappedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given timeDependentAverageMappedFvPatchField
        //  onto a new patch
        timeDependentAverageMappedFvPatchField
        (
            const timeDependentAverageMappedFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        timeDependentAverageMappedFvPatchField
        (
            const timeDependentAverageMappedFvPatchField<Type>&
        ) = delete;

        //- Copy constructor setting internal field reference
        timeDependentAverageMappedFvPatchField
        (
            const timeDependentAverageMappedFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new timeDependentAverageMappedFvPatchField<Type>(*this, iF)
            );
        }


    // Member Functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            //  Used to update fields following mesh topology change
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            //  Used to reconstruct fields
            virtual void rmap(const fvPatchField<Type>&, const labelList&);

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "timeDependentAverageMappedFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
