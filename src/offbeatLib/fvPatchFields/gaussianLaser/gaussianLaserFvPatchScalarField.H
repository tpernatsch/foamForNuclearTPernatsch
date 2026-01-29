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
    Foam::gaussianLaserFvPatchScalarField

Description
    BOundary condition to prescribe a gaussian LASER intensity profile on a 
    patch by providing LASER power, beam waist and coefficient of reflection.

Usage
    fuelTop
    {
        type            gaussianLaser;

        // Uncomment this for const laser power
        // power           150;

        // Time dependent laser power
        powerList
        {
            values
            (
                // (time in (s) - power in (W))
                (0.0    0) 
                (60     150)
                (1e9    150)
            );

            outOfBounds         clamp; // optional out-of-bounds handling
            interpolationScheme linear; // optional interpolation method
        }

        beamWaist       3.5e-3; // beam waist in m
        reflectionCoeff 0.20; // reflection coeff (1-emissivity)
        value           $internalField;
    }
 
SourceFiles
    gaussianLaserFvPatchScalarField.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
    Switzerland, Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro, M.Reymond - EPFL\n

\date 
    May 2024
 
\*---------------------------------------------------------------------------*/
 
#ifndef gaussianLaserFvPatchScalarField_H
#define gaussianLaserFvPatchScalarField_H
 
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
        Class gaussianLaserFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/
 
class gaussianLaserFvPatchScalarField
:
    public fixedValueFvPatchField<scalar>
{
    // Private Data

protected:

    // Protected Data

        //- LASER power (W)
        mutable scalar P_;

        //- LASER beam waist
        const scalar w_;

        //- Coefficient of reflection
        const scalar R_;

#ifdef OPENFOAMFOUNDATION
        //- Time varying list of LASER powers (W)
        autoPtr<Function1s::Table<scalar> > powerList_;
#elif OPENFOAMESI
        typedef Foam::Function1Types::Table<scalar> scalarTable;

        //- Time varying list of LASER powers (W)
        autoPtr<scalarTable> powerList_;
#endif        

    // Protected Member Functions
 
public:
 
    //- Runtime type information
    TypeName("gaussianLaser");
 
 
    // Constructors
 
        //- Construct from patch and internal field
        gaussianLaserFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );
 
        //- Construct from patch, internal field and dictionary
        gaussianLaserFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );
 
        //- Construct by mapping the given gaussianLaserFvPatchScalarField
        //  onto a new patch
        gaussianLaserFvPatchScalarField
        (
            const gaussianLaserFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&,
            const bool mappingRequired=true
        );
 
 
        //- Copy constructor setting internal field reference
        gaussianLaserFvPatchScalarField
        (
            const gaussianLaserFvPatchScalarField&,
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
                new gaussianLaserFvPatchScalarField(*this, iF)
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