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
     Foam::fixedDisplacementFvPatchVectorField
 
Description
    fixedDisplacement patchfield. 
    The fixed displacement can be provided as a time series.

Usage
    \verbatim
    cladOuter
    {
        type            fixedDisplacement;
        displacementSeries
        {
            type    table;

            values
            (
                (0.0 (0 0 0))
                (100.0 (0.1 0 0.1))
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
        type            fixedDisplacement;
        value           uniform (0.1 0 0.1);
    }    
    \endverbatim   
 
SourceFiles
    fixedDisplacementFvPatchVectorField.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
    Switzerland, Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021
 
\*---------------------------------------------------------------------------*/
 
#ifndef fixedDisplacementFvPatchVectorField_H
#define fixedDisplacementFvPatchVectorField_H
 
#include "fixedValueFvPatchField.H"
#include "axialProfile.H"
#include "InterpolateTables.H"

#include "primitiveFields.H"
#include "objectRegistry.H"
#include "fvMesh.H"
#include "volFields.H"
#include "Table.H"

#ifdef OPENFOAMESI
    #include "Function1.H"
#endif   
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
namespace Foam
{
 
/*---------------------------------------------------------------------------*\
        Class fixedDisplacementFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/
 
class fixedDisplacementFvPatchVectorField
:
    public fixedValueFvPatchField<vector>
{
    // Private Data

protected:

    // Protected Data

        //- Total displacement field (D not DD)
        vectorField totalDisp_;

#ifdef OPENFOAMFOUNDATION
        //- Displacement time series
        autoPtr<Function1s::Table<vector> > dispSeries_;
#elif OPENFOAMESI
        //- Displacement time series
        autoPtr<Function1<vector> >  dispSeries_;
#endif 

    // Protected Member Functions
 
public:
 
    //- Runtime type information
    TypeName("fixedDisplacement");
 
 
    // Constructors
 
        //- Construct from patch and internal field
        fixedDisplacementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );
 
        //- Construct from patch, internal field and dictionary
        fixedDisplacementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&,
            const bool valueRequired=true
        );
 
        //- Construct by mapping the given fixedDisplacementFvPatchVectorField
        //  onto a new patch
        fixedDisplacementFvPatchVectorField
        (
            const fixedDisplacementFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&,
            const bool mappingRequired=true
        );
 
 
        //- Copy constructor setting internal field reference
        fixedDisplacementFvPatchVectorField
        (
            const fixedDisplacementFvPatchVectorField&,
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
                new fixedDisplacementFvPatchVectorField(*this, iF)
            );
        }
 
 
    // Member Functions

        //- Write
        virtual void write(Ostream&) const;
    
    // Evaluation functions
 
        //- Return the matrix diagonal coefficients corresponding to the
        //  evaluation of the value of this patchField with given weights
        virtual void updateCoeffs ();

 };
 
 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
 } // End namespace Foam
 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
 #endif
 
 // ************************************************************************* //