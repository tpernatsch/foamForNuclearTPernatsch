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
    Foam::timeDependentAxialAzimuthalProfileFvPatchScalarField
 
Description
    Boundary condition to specify a time-dependent axial and azimuthal temperature profile.
    Assumes the axial and azimuthal profiles are independent, i.e.
        T(z,theta,t) = f_z*(z,t) + f_theta(theta,t)

    \htmlinclude fvPatchField_timeDependentAxialAzimuthalProfileFvPatchScalarField.html
 
SourceFiles
    timeDependentAxialAzimuthalProfileFvPatchScalarField.C

\todo

\mainauthor
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\contribution

\date 
    January 2024
 
\*---------------------------------------------------------------------------*/
 
#ifndef timeDependentAxialAzimuthalProfileFvPatchScalarField_H
#define timeDependentAxialAzimuthalProfileFvPatchScalarField_H
 
#include "timeDependentAxialProfileFvPatchScalarField.H"
#include "InterpolateTables.H"
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
namespace Foam
{
 
/*---------------------------------------------------------------------------*\
        Class timeDependentAxialAzimuthalProfileFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/
 
class timeDependentAxialAzimuthalProfileFvPatchScalarField
:
    public timeDependentAxialProfileFvPatchScalarField
{
    // Private Data

    // Private Member Functions
        
        //- Calculate the angle at each boundary face
        void calculateTheta();
    
protected:

    // Protected Data
    
		//- Tabulates theta values
		scalarField thetaValues_;
		
		//- Azimuthal tabulated data_;
		FieldField<Field, scalar> thetaData_;
		
        //- Axial location interpolation method
        interpolateTableBase::interpolationMethod thetaMethod_;
        
        //- Interpolation table
        scalarFieldInterpolateTable thetaTable_;
        
		//- Radial-axis from which to calculate theta
		vector rAxis_;
        
        //- Theta values on boundary faces
        scalarField theta_;

public:
 
    //- Runtime type information
    TypeName("timeDependentAxialAzimuthalProfileT");
 
 
    // Constructors
 
        //- Construct from patch and internal field
        timeDependentAxialAzimuthalProfileFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );
 
        //- Construct from patch, internal field and dictionary
        timeDependentAxialAzimuthalProfileFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&,
            const bool valueRequired=true
        );
 
        //- Construct by mapping the given timeDependentAxialAzimuthalProfileFvPatchScalarField
        //  onto a new patch
        timeDependentAxialAzimuthalProfileFvPatchScalarField
        (
            const timeDependentAxialAzimuthalProfileFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );
 
 
        //- Copy constructor setting internal field reference
        timeDependentAxialAzimuthalProfileFvPatchScalarField
        (
            const timeDependentAxialAzimuthalProfileFvPatchScalarField&,
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
                new timeDependentAxialAzimuthalProfileFvPatchScalarField(*this, iF)
            );
        }
 
    // Member Functions

    // Evaluation functions
 
        //- Return the matrix diagonal coefficients corresponding to the
        //  evaluation of the value of this patchField with given weights
        virtual void updateCoeffs();

 };
 
 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
 } // End namespace Foam
 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
 #endif
 
 // ************************************************************************* //
