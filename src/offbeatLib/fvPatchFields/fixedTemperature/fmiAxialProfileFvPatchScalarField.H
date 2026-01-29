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
    Foam::fmiAxialProfileFvPatchScalarField
 
Description
    Boundary condition to impose a fixed axial temperature profile from an FMI.
    If the oxideLayer is present, an additional thermal resistance corresponding
    to the oxide is considered in order to set the patch temperature.

Usage
    \verbatim
    claddingOuterSurface
    {
        type            fmiAxialProfileT;
        axialProfileDict
        {
            axialLocationsNameFromFMU   offbeatAxialLocations;
            axialProfileNameFromFMU     offbeatTstructAxial;
            axialInterpolationMethod    linear;
            initialData                 7{300};
            // NOTE : axial locations are NOT normalized
            initialAxialLocations       ( 0.0 0.5 1.0 1.5 2.0 2.5 3.0 );
        }
        value           uniform 300;
    }
    // or
    claddingOuterSurface // Imposed a uniform axial distribution
    {
        type            fmiAxialProfileT;
        axialProfileDict
        {
            axialLocationsNameFromFMU   offbeatAxialLocations;
            axialProfileNameFromFMU     offbeatTstructAxial;
            axialInterpolationMethod    linear;
            initialData                 300;
        }
        value           uniform 300;
    }
    \endverbatim
 
SourceFiles
    fmiAxialProfileFvPatchScalarField.C

\mainauthor
    Thomas Guilbaud, EPFL/Transmutex SA, Switzerland

\date 
    October 2023
 
\*---------------------------------------------------------------------------*/
 
#ifndef fmiAxialProfileFvPatchScalarField_H
#define fmiAxialProfileFvPatchScalarField_H
 
#include "fixedTemperatureFvPatchScalarField.H"
#include "axialProfile.H"
#include "InterpolateTables.H"

#include "globalOptions.H"
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
namespace Foam
{
 
/*---------------------------------------------------------------------------*\
        Class fmiAxialProfileFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/
 
class fmiAxialProfileFvPatchScalarField
:
    public fixedTemperatureFvPatchScalarField
{
    // Private Data

        //- FMI port name for axial locations
        word axialLocationsNameFromFMU_;
    
        //- FMI port name for axial profile
        word axialProfileNameFromFMU_;
    
        //- Axial location interpolation method
        interpolateTableBase::interpolationMethod zMethod_;

        //- Reference axial locations
        scalarField zRef_;
        
        //- Axial profile dict
        dictionary dict_;
 
public:
 
    //- Runtime type information
    TypeName("fmiAxialProfileT");
 
 
    // Constructors
 
        //- Construct from patch and internal field
        fmiAxialProfileFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );
 
        //- Construct from patch, internal field and dictionary
        fmiAxialProfileFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&,
            const bool valueRequired=true
        );
 
        //- Construct by mapping the given fmiAxialProfileFvPatchScalarField
        //  onto a new patch
        fmiAxialProfileFvPatchScalarField
        (
            const fmiAxialProfileFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );
 
 
        //- Copy constructor setting internal field reference
        fmiAxialProfileFvPatchScalarField
        (
            const fmiAxialProfileFvPatchScalarField&,
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
                new fmiAxialProfileFvPatchScalarField(*this, iF)
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