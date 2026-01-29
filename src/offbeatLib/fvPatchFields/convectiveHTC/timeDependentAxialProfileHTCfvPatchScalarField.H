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
    The `timeDependentAxialProfileHTC` boundary condition allows to simulate an heat
    sink boundary defined by a convective heat transfer coefficient \f$h\f$ and a
    sink temperature \f$T_0\f$ that vary axially and during time.
    
    At each iteration the temperature at the boundary is calculated using the classic
    convective heat transfer formula as done in the parent class `convectiveHTC`.
    The only difference is that the both the sink temperature and the heat trasnsfer
    coefficient can vary axially and in time.

    \note
    The axial locations provided to this boundary conditions are NOT normalized, thus they
    are absolute values expressed in \f$m\f$.

________________________________________________________________________________

\par Options
    The `timeDependentAxialProfileHTC` fvPatchField can be selected in the patch 
    subdictionary inside the `boundaryField` subdictionary of the temperature field.

    Parameters in the patch subdictionary for `timeDependentAxialProfileHTC`:
    
    - <b>`axialProfileDict`</b> - A sub-dictionary that contains the data for the 
    \f$T_0\f$ and \f$h\f$ axial- and time-dependent profiles. See following parameters
    section.
    - <b>`kappa`</b> - The name of the conductivity field. **By default, it is set
    to "k", which is the field name for conductivity used by the OFFBEAT material class.**
    - <b>`value`</b> - Specifies the initial value of the patch temperature field.

    Parameters in the `axialProfileDict` subdictionary:

    - <b>`timePoints`</b> - A list of time values at which the axial profiles are provided. 
    The time unit depend on the userTime selected by the user (seconds by default).
    - <b>`axialLocations`</b> - A list of axial locations (in `m`) at which the sink temperature
    and heat exchange coefficients are provided (for each time point). **Axial locations
    are not normalized!**
    - <b>`hData`</b> - 2D table for the heat trasnfer coefficient axial profiles (in W/m\f$^2\f$/K).
    Axial locations are on the x axis, time is on the y axis (i.e. every row provides the axial 
    profile at a given time point in the `timePoints` list).
    - <b>`T0Data`</b> - 2D table for the sink temperature axial profiles (in K) . Axial
    locations are on the x axis, time is on the y axis (i.e. every row provides the axial profile
    at a given time point in the `timePoints` list).
    - <b>`timeInterpolationMethod`</b> - It selects the time 
    interpolation method for time steps that fall in between the time points indicated in 
    the `timePoints` list.**Either `step` or `linear` interpolation can be selected.**
    - <b>`axialInterpolationMethod`</b> - It selects the axial 
    interpolation method for the axial profile. **Either `step` or `linear` interpolation
    can be selected.**

________________________________________________________________________________

Usage
    To apply the `timeDependentAxialProfileHTC` boundary condition to a given patch,
    the following example can be used as a template:
    
    \verbatim
        cladOuter
        {
            type            timeDependentAxialProfileHTC;

            axialProfileDict
            {
                timePoints     (0 1e3 1.261e8);
                axialLocations (0.0 1.0 2.0 3.2);
                hData          (
                    (5e3 5.01e3 5.02e3 5.01e3)
                    (5e3 5.01e3 5.02e3 5.01e3)
                    (5e3 5.01e3 5.02e3 5.01e3)
                );
                T0Data          (
                    (350 370 380 365)
                    (360 380 370 360)
                    (340 360 390 376)
                );

                axialInterpolationMethod "linear";
                timeInterpolationMethod "linear";
            }

            value           $internalField;
        }
    \endverbatim
    
    <!-- </details> -->
        

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
