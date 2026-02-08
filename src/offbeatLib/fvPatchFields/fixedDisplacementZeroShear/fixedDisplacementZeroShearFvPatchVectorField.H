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
    Foam::fixedDisplacementZeroShearFvPatchVectorField

Description
    The `fixedDisplacementZeroShear` boundary condition fixes the displacement 
    component normal to the patch, while the shear traction is set to zero. 
    This fvPatchField is similar to a symmetry boundary condition but allows 
    the material to slide parallel to the patch.

    Optionally, the fixed normal displacement can be specified as a time series.

    \note
    This fvPatchField is taken from [solids4Foam](https://www.solids4foam.com/) library.

________________________________________________________________________________

\par Options
    The `fixedDisplacementZeroShear` fvPatchField can be selected in the patch subdictionary 
    inside the `boundaryField` subdictionary of the displacement field.

    Parameters in the patch subdictionary for `fixedDisplacementZeroShear`:

    - <b>`displacementSeries`</b> - A subdictionary specifying the time series of 
    displacement values.
    - <b>`value`</b> - Only used if `displacementSeries` is not present, it sets the
    fixed displacement value (constant in time).

    Parameters in the subdictionary `displacementSeries`:
    
    - <b>`type`</b> - Specifies the type of the series, e.g., `table`.
    - <b>`values`</b> - A list of time-value pairs defining the displacement series. 
    The displacement is specified in vector form for each time point.
    - <b>`outOfBounds`</b> - Optional keyword to handle values outside the defined
    time range. Possible values: `clamp` or `extrapolate`. **Default: `clamp`.**
    - <b>`interpolationScheme`</b> - Optional keyword to define the interpolation
    method between time points. Possible values: `linear` or `step`. **Default: `linear`.**

    \warning
    Displacement values must always be provided as vectors, but only the normal
    component (normal to the patch) will be applied.

________________________________________________________________________________

Usage
    To apply the `fixedDisplacementZeroShear` boundary condition to a given patch, 
    the following example can be used as a template:

    \verbatim
        cladOuter
        {
            type            fixedDisplacementZeroShear;

            // For time-dependent series
            displacementSeries
            {
                type    table;

                values
                (
                    (0.0 (0 0 0))
                    (100.0 (0 0 0.1))
                );

                outOfBounds         clamp;     // optional out-of-bounds handling
                interpolationScheme linear;    // optional interpolation method
            }

            // // Alternatively, for a constant normal displacement
            // value           uniform (0 0 0.1);
        }    
    \endverbatim       

SourceFiles
    fixedDisplacementZeroShearFvPatchVectorField.C

\mainauthor
    Alessandro Scolaro - EFPL (derived from the work of Philip Cardiff, UCD/UT)


\*---------------------------------------------------------------------------*/

#ifndef fixedDisplacementZeroShearFvPatchVectorField_H
#define fixedDisplacementZeroShearFvPatchVectorField_H

#include "fvPatchFields.H"
#include "solidDirectionMixedFvPatchVectorField.H"
#include "Switch.H"
#include "Table.H"

#ifdef OPENFOAMESI
    #include "Function1.H"
#endif   

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class fixedDisplacementZeroShearFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class fixedDisplacementZeroShearFvPatchVectorField
:
    public solidDirectionMixedFvPatchVectorField
{
protected:
    
    // Protected data

        //- Constant total displacement field
        vectorField totalDisp_;
        
#ifdef OPENFOAMFOUNDATION
        //- Displacement time series
        autoPtr<Function1s::Table<vector>> dispSeries_;
#elif OPENFOAMESI
        //- Displacement time series
        autoPtr<Function1<vector> >  dispSeries_;
#endif        


public:

    //- Runtime type information
    TypeName("fixedDisplacementZeroShear");


    // Constructors

        //- Construct from patch and internal field
        fixedDisplacementZeroShearFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedDisplacementZeroShearFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&,
            const bool valueRequired = true
        );

        //- Construct by mapping given
        // fixedDisplacementZeroShearFvPatchVectorField
        //  onto a new patch
        fixedDisplacementZeroShearFvPatchVectorField
        (
            const fixedDisplacementZeroShearFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&,
            const bool mappingRequired=true
        );

        //- Construct as copy setting internal field reference
        fixedDisplacementZeroShearFvPatchVectorField
        (
            const fixedDisplacementZeroShearFvPatchVectorField&,
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
                new fixedDisplacementZeroShearFvPatchVectorField(*this, iF)
            );
        }

    // Destructor

        virtual ~fixedDisplacementZeroShearFvPatchVectorField()
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
