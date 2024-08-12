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
    Foam::xzTemperatureProfileFvPatchScalarField
 
Description

    Boundary condition to apply a temperature profile constant in time but 
    varying along x and z directions.

    \htmlinclude fvPatchField_xzTemperatureProfileFvPatchScalarField.html


SourceFiles
    xzTemperatureProfileFvPatchScalarField.C

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
 
#ifndef xzTemperatureProfileFvPatchScalarField_H
#define xzTemperatureProfileFvPatchScalarField_H
 
#include "fixedTemperatureFvPatchScalarField.H"
#include "axialProfile.H"
#include "InterpolateTables.H"

#include "globalOptions.H"
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
namespace Foam
{
 
/*---------------------------------------------------------------------------*\
        Class xzTemperatureProfileFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/
 
class xzTemperatureProfileFvPatchScalarField
:
    public fixedTemperatureFvPatchScalarField
{
    // Private Data
    
        //- Tabulated z locations
        scalarField zValues_;

        //- Tabulated x locations
        scalarField xValues_;
        
        //- Tabulated Profile data
        FieldField<Field, scalar> profileData_;
    
        //- z location interpolation method
        interpolateTableBase::interpolationMethod zMethod_;
        
        //- x location interpolation method
        interpolateTableBase::interpolationMethod xMethod_;
        
        //- Interpolation table
        scalarFieldInterpolateTable table_;
        
        //- Axial profile dict
        dictionary dict_;
 
public:
 
    //- Runtime type information
    TypeName("xzTemperatureProfile");
 
 
    // Constructors
 
        //- Construct from patch and internal field
        xzTemperatureProfileFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );
 
        //- Construct from patch, internal field and dictionary
        xzTemperatureProfileFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&,
            const bool valueRequired=true
        );
 
        //- Construct by mapping the given xzTemperatureProfileFvPatchScalarField
        //  onto a new patch
        xzTemperatureProfileFvPatchScalarField
        (
            const xzTemperatureProfileFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&,
            const bool mappingRequired=true
        );
 
 
        //- Copy constructor setting internal field reference
        xzTemperatureProfileFvPatchScalarField
        (
            const xzTemperatureProfileFvPatchScalarField&,
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
                new xzTemperatureProfileFvPatchScalarField(*this, iF)
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