/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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
    Foam::convectiveHTCFvPatchScalarField

Description
    The `convectiveHTC` boundary condition allows to simulate an heat sink boundary
    defined by a convective heat transfer coefficient \f$h\f$ and a sink temperature
    \f$T_0\f$. 
    
    This fvPatchField computes the temperature \f$T_b\f$ to be applied to the
    patch using the following expression for heat flux \f$q"\f$:
    
    $$
    \begin{aligned}
    q" = h\cdot(T_b-T_0)
    \end{aligned}
    $$
    
    where \f$T_b\f$ is the boundary temperature and \f$q"\f$ is computed from 
    the Fourier's Law as:
    
    $$
    \begin{aligned}
    q" = -k\cdot \nabla T
    \end{aligned}
    $$
    
    with \f$k\f$ as the thermal conductivity.

________________________________________________________________________________

\par Options
    The `convectiveHTC` fvPatchField can be selected in the patch subdictionary
    inside the `boundaryField` subdictionary of the temperature field.

    Parameters in the patch subdictionary for `convectiveHTC`:

    - <b>`T0`</b> - Fixed sink temperature in K.
    - <b>`h`</b> - Fixed heat conductance in W/m\f$^2\f$/K.
    - <b>`kappa`</b> - The name of the conductivity field. **By default, it is set
    to "k", which is the field name for conductivity used by the OFFBEAT material class.**
    - <b>`value`</b> - Specifies the initial value of the patch temperature field.
    
________________________________________________________________________________

Usage
    To apply the `convectiveHTC` boundary condition to a given patch, the following
    example can be used as a template:
    
    \verbatim
        cladOuter
        {
            type            convectiveHTC;

            // Sink temperature in K
            T0              uniform 600;

            // Heat transfer coefficient in W/m^2/K
            h               uniform 50000;

            value           uniform 300;
        }
    \endverbatim	    

SourceFiles
    convectiveHTCFvPatchScalarField.C

\mainauthor
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)    

\contribution
    A. Scolaro, E. Brunetto, C. Fiorina - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
    Switzerland, Laboratory for Reactor Physics and Systems Behaviour)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef convectiveHTCFvPatchScalarField_H
#define convectiveHTCFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class convectiveHTCFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class convectiveHTCFvPatchScalarField
:
    public fvPatchField<scalar>
{
    // Private data

        //- Convective heat transfer coefficient
        scalarField h_;
        
        //- Environment temperature
        scalarField T0_;
        
        //- Name of thermal conductivity field
        word kappaName_;
        
        //- Boundary heat transfer factor - working variable
        scalarField alpha_;
        
public:

    //- Runtime type information
    TypeName("convectiveHTC");


    // Constructors

        //- Construct from patch and internal field
        convectiveHTCFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        convectiveHTCFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given convectiveHTCFvPatchScalarField
        //  onto a new patch
        convectiveHTCFvPatchScalarField
        (
            const convectiveHTCFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        convectiveHTCFvPatchScalarField
        (
            const convectiveHTCFvPatchScalarField&,
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
                new convectiveHTCFvPatchScalarField(*this, iF)
            );
        }


    // Member functions
    
        //- Access the heat transfer coefficient
        inline scalarField& h()
        {
            return h_;
        }
        
        //- Access the heat transfer coefficient
        inline const scalarField& h() const
        {
            return h_;
        }
        
        //- Access the environment temperature
        inline scalarField& T0()
        {
            return T0_;
        }

        //- Access the environment temperature
        inline const scalarField& T0() const
        {
            return T0_;
        }

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
