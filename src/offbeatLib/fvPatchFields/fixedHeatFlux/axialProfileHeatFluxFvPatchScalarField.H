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
    Foam::axialProfileHeatFluxFvPatchScalarField

Description
    Fixed heat flux boundary condition.

    \htmlinclude fvPatchField_axialProfileHeatFluxFvPatchScalarField.html

SourceFiles
    axialProfileHeatFluxFvPatchScalarField.C

\mainauthor
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)    

\date 
    November 2023

\*---------------------------------------------------------------------------*/

#ifndef axialProfileHeatFluxFvPatchScalarField_H
#define axialProfileHeatFluxFvPatchScalarField_H

#include "fixedHeatFluxFvPatchScalarField.H"
#include "axialProfile.H"
#include "InterpolateTables.H"
#include "globalOptions.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class axialProfileHeatFluxFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class axialProfileHeatFluxFvPatchScalarField
:
    public fixedHeatFluxFvPatchScalarField
{
    // Private data

        typedef interpolateTableBase::interpolationMethod interpolationMethod;
    
        //- Axial profile dict
        dictionary profileDict_;
		
		//- Reference mesh
		const fvPatch& referencePatch_;

        //- List of time points
        scalarField timePoints_;

        //- List of time points
        scalarField axialPoints_;
        
        //- List of average values
        scalarField averageData_;

        //- List of average values
        FieldField<Field, scalar> axialData_;

        //- Time interpolation method
        interpolateTableBase::interpolationMethod tMethod_;
        
        //- Axial location interpolation method
        interpolateTableBase::interpolationMethod zMethod_;
        
        //- Time-dependent table for average heat flux
        scalarInterpolateTable averageTable_;

        //- Time-dependent table for axial profile
        scalarFieldInterpolateTable axialTable;

		//- Z-axis
		vector axis_;

        //- Axial locations of boundary faces
        scalarField z_;
        
    // Private Functions
    
        //- Calculate the axial location of each boundary face
        void calculateAxialLocations();
        
public:

    //- Runtime type information
    TypeName("axialProfileHeatFlux");


    // Constructors

        //- Construct from patch and internal field
        axialProfileHeatFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        axialProfileHeatFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given axialProfileHeatFluxFvPatchScalarField
        //  onto a new patch
        axialProfileHeatFluxFvPatchScalarField
        (
            const axialProfileHeatFluxFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        axialProfileHeatFluxFvPatchScalarField
        (
            const axialProfileHeatFluxFvPatchScalarField&,
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
                new axialProfileHeatFluxFvPatchScalarField(*this, iF)
            );
        }


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
                const fvPatchScalarField&,
                const labelList&
            );
        
            
        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write to an output stream
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
