/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::velocityRundownFvPatchVectorField

Group
    grpInletBoundaryConditions

Description
    Velocity inlet boundary condition either correcting the extrapolated
    velocity or creating a uniform velocity field normal to the patch adjusted
    to match the specified flow rate

    For a mass-based flux:
    - the flow rate should be provided in kg/s
    - if \c rho is "none" the flow rate is in m3/s
    - otherwise \c rho should correspond to the name of the density field
    - if the density field cannot be found in the database, the user must
      specify the inlet density using the \c rhoInlet entry

    For a volumetric-based flux:
    - the flow rate is in m3/s

Usage
    Sets a time-dependent velocity in the form U(t) = U0*((A+B*(t-t0))^C)

See also
    Foam::fixedValueFvPatchField

SourceFiles
    velocityRundownFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef velocityRundownFvPatchVectorField_H
#define velocityRundownFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class velocityRundownFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class velocityRundownFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        vectorField U0_;
        scalar A_;
        scalar B_;
        scalar C_;
        scalar t0_;

public:

   //- Runtime type information
   TypeName("velocityRundown");


   // Constructors

        //- Construct from patch and internal field
        velocityRundownFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        velocityRundownFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  velocityRundownFvPatchVectorField
        //  onto a new patch
        velocityRundownFvPatchVectorField
        (
            const velocityRundownFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        velocityRundownFvPatchVectorField
        (
            const velocityRundownFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new velocityRundownFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        velocityRundownFvPatchVectorField
        (
            const velocityRundownFvPatchVectorField&,
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
                new velocityRundownFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
