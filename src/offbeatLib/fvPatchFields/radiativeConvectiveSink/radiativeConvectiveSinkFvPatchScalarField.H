/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::radiativeConvectiveSinkFvPatchScalarField

Description
    Heat sink boundary for solid temperature defined by radiative and
    convective heat transfer coefficients.

________________________________________________________________________________

\par Options
    The `radiativeConvectiveSink` fvPatchField can be selected in the patch subdictionary
    inside the `boundaryField` subdictionary of the temperature field.

    Parameters in the patch subdictionary for `radiativeConvectiveSink`:

    - <b>`kappa`</b> - The name of the conductivity field. **By default, it is set
    to "k", which is the field name for conductivity used by the OFFBEAT material class.**    
    - <b>`alpha`</b> - Specifies the value of the gap conductance field in \f$W/K/m^2\f$.
    This can either be a `nonUniformList<scalar>` with a value defined for each face,
    or a `uniform` field with a single value applied across all faces.    
    - <b>`emissivity`</b> - Specifies the value of the emissivity field.
    This can either be a `nonUniformList<scalar>` with a value defined for each face,
    or a `uniform` field with a single value applied across all faces.
    - <b>`value`</b> - Specifies the initial value of the patch temperature field.

________________________________________________________________________________

Usage
    \verbatim
        cladOuter
        {
            type            radiativeConvectiveSink;

            // Sink temperature
            T0              uniform 600;
            // Convective heat transfer coefficient
            alpha           uniform 50000;
            // Radiative emissivity
            emissivity      uniform 0.5;

            value           $internalField;
        }       
    \endverbatim    

SourceFiles
    radiativeConvectiveSinkFvPatchScalarField.C

\mainauthor
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)
    

\contribution
    A. Scolaro, E. Brunetto, C. Fiorina - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
    Switzerland, Laboratory for Reactor Physics and Systems Behaviour)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef radiativeConvectiveSinkFvPatchScalarField_H
#define radiativeConvectiveSinkFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class radiativeConvectiveSinkFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class radiativeConvectiveSinkFvPatchScalarField
:
    public fvPatchField<scalar>
{
    // Private data

        //- Convective heat transfer coefficient
        scalarField H_;
        
        //- Radiative emissivity
        scalarField emissivity_;
        
        //- Sink temperature
        scalarField T0_;
        
        //- Name of thermal conductivity field
        word kappaName_;
        
        //- Boundary heat transfer factor - working variable
        scalarField alpha_;
        
public:

    //- Runtime type information
    TypeName("radiativeConvectiveSink");


    // Constructors

        //- Construct from patch and internal field
        radiativeConvectiveSinkFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        radiativeConvectiveSinkFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given radiativeConvectiveSinkFvPatchScalarField
        //  onto a new patch
        radiativeConvectiveSinkFvPatchScalarField
        (
            const radiativeConvectiveSinkFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&,
            const bool mappingRequired=true
        );

        //- Construct as copy setting internal field reference
        radiativeConvectiveSinkFvPatchScalarField
        (
            const radiativeConvectiveSinkFvPatchScalarField&,
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
                new radiativeConvectiveSinkFvPatchScalarField(*this, iF)
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
