/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::zeroCurrentActinidesRedistributionFvPatchScalarField

Description


Usage

SourceFiles
    zeroCurrentActinidesRedistributionFvPatchScalarField.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
    Switzerland, Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro, C. Fiorina - EPFL

\date 
    October 2022

\*---------------------------------------------------------------------------*/

#ifndef zeroCurrentActinidesRedistributionFvPatchScalarField_H
#define zeroCurrentActinidesRedistributionFvPatchScalarField_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class zeroCurrentActinidesRedistributionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class zeroCurrentActinidesRedistributionFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{

    // Private Data
        
protected:

    // Protected Data
  
        //- Universal gas constant [J/K/mol)]
        const scalar R_;
        
        //- Molar heat of transport [J/mol]
        const scalar Q_;
        
        //- Pores diameter [m]
        const scalar poreDiameter_;

        //- Pores thickness  [m]
        const scalar poreThickness_;

        //- Coefficient in pore migration term [1/K]
        const scalar A_;
        
public:

    //- Runtime type information
    TypeName("zeroCurrentActinidesRedistribution");


    // Constructors

        //- Construct from patch and internal field
        zeroCurrentActinidesRedistributionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        zeroCurrentActinidesRedistributionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  zeroCurrentActinidesDiffusion onto a new patch
        zeroCurrentActinidesRedistributionFvPatchScalarField
        (
            const zeroCurrentActinidesRedistributionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        zeroCurrentActinidesRedistributionFvPatchScalarField
        (
            const zeroCurrentActinidesRedistributionFvPatchScalarField&,
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
                new zeroCurrentActinidesRedistributionFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }



    // Member functions

        // Access

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap(const fvPatchScalarField&, const labelList&);


        //- Set the gradient
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
