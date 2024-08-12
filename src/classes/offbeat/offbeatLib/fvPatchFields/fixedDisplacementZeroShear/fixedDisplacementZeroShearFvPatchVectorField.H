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
    The component of the displacement normal to the patch is applied,
    and the shear traction is set to zero.

    Patch contains non-orthogonal correction.

    Optionally time-varying by specifying a time series.

Usage
    \verbatim
    cladOuter
    {
        type            fixedDisplacementZeroShear;
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
    }    
    \endverbatim  

    or with constant fixed value:

    \verbatim
    cladOuter
    {
        type            fixedDisplacementZeroShear;
        value           uniform (0 0 0.1);
    }    
    \endverbatim         

SourceFiles
    fixedDisplacementZeroShearFvPatchVectorField.C

\mainauthor
    Philip Cardiff, UCD/UT. All rights reserved. Based on velocityDisplacement by Aleksandar Karac


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
            const fvPatchFieldMapper&
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
