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
    Foam::timeDependentAxialProfileFvPatchScalarField
 
Description
    Boundary condition to impose a time-dependent axial temperature profile.

    \htmlinclude fvPatchField_timeDependentAxialProfileFvPatchScalarField.html
 
SourceFiles
    timeDependentAxialProfileFvPatchScalarField.C

\todo
    current version of write(os) is not ideal since it writes the entire 
    lists for zValues, timeValues and profileData in the boundaryField section 
    of T

\mainauthor
    E. Brunetto, A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
    Switzerland, Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021
 
\*---------------------------------------------------------------------------*/
 
#ifndef timeDependentAxialProfileFvPatchScalarField_H
#define timeDependentAxialProfileFvPatchScalarField_H
 
#include "fixedTemperatureFvPatchScalarField.H"
#include "axialProfile.H"
#include "InterpolateTables.H"

#include "globalOptions.H"
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
namespace Foam
{
 
/*---------------------------------------------------------------------------*\
        Class timeDependentAxialProfileFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/
 
class timeDependentAxialProfileFvPatchScalarField
:
    public fixedTemperatureFvPatchScalarField
{
    // Private Data
    
        //- Tabulated axial locations
        scalarField zValues_;
        
        //- Tabulated time values
        const scalarField timeValues_;
        
        //- Tabulated Profile data
        FieldField<Field, scalar> profileData_;
    
        //- Axial location interpolation method
        interpolateTableBase::interpolationMethod zMethod_;
        
        //- Time interpolation method
        interpolateTableBase::interpolationMethod tMethod_;
        
        //- Interpolation table
        scalarFieldInterpolateTable table_;

        //- Reference axial locations
        scalarField zRef_;
        
        //- Axial profile dict
        dictionary dict_;
 
public:
 
    //- Runtime type information
    TypeName("timeDependentAxialProfileT");
 
 
    // Constructors
 
        //- Construct from patch and internal field
        timeDependentAxialProfileFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );
 
        //- Construct from patch, internal field and dictionary
        timeDependentAxialProfileFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&,
            const bool valueRequired=true
        );
 
        //- Construct by mapping the given timeDependentAxialProfileFvPatchScalarField
        //  onto a new patch
        timeDependentAxialProfileFvPatchScalarField
        (
            const timeDependentAxialProfileFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&,
            const bool mappingRequired=true
        );
 
 
        //- Copy constructor setting internal field reference
        timeDependentAxialProfileFvPatchScalarField
        (
            const timeDependentAxialProfileFvPatchScalarField&,
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
                new timeDependentAxialProfileFvPatchScalarField(*this, iF)
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