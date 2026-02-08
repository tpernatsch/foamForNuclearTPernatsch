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
    Foam::fixedHeatFluxFvPatchScalarField

Description
    Fixed heat flux boundary condition.

    \htmlinclude fvPatchField_fixedHeatFluxFvPatchScalarField.html

SourceFiles
    fixedHeatFluxFvPatchScalarField.C

\mainauthor
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)    

\date 
    November 2023

\*---------------------------------------------------------------------------*/

#ifndef fixedHeatFluxFvPatchScalarField_H
#define fixedHeatFluxFvPatchScalarField_H

#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class fixedHeatFluxFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class fixedHeatFluxFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{

protected:

    // Protected data

        //- Heat Flux
        scalarField q_;
        
        //- Name of thermal conductivity field
        word kappaName_;
        
public:

    //- Runtime type information
    TypeName("fixedHeatFlux");


    // Constructors

        //- Construct from patch and internal field
        fixedHeatFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedHeatFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedHeatFluxFvPatchScalarField
        //  onto a new patch
        fixedHeatFluxFvPatchScalarField
        (
            const fixedHeatFluxFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        fixedHeatFluxFvPatchScalarField
        (
            const fixedHeatFluxFvPatchScalarField&,
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
                new fixedHeatFluxFvPatchScalarField(*this, iF)
            );
        }


    // Member functions
    
        //- Access the heat transfer coefficient
        inline scalarField& q()
        {
            return q_;
        }
        
        //- Access the heat transfer coefficient
        inline const scalarField& q() const
        {
            return q_;
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

        //- Write to an output stream
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
