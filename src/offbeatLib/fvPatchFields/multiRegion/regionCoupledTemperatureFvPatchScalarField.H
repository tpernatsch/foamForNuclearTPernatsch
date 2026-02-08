/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::regionCoupledTemperatureFvPatchScalarField

Description
    turbulentTemperatureCoupledBaffleMixedFvPatchScalarField for OFFBEAT
    usage.
    Mixed boundary condition for temperature, to be used for heat-transfer
    with another region in a CHT case.

    Specifies gradient and temperature such that the equations are the same
    on both sides:
      - refGradient = qs_/kappa
      - refValue = neighbour value
      - valueFraction = kappaByDeltaNbr/(kappaByDeltaNbr + kappaByDelta)

    where kappaByDelta is heat-transfer coefficient kappa*deltaCoeffs.

Usage
    \table
        Property     | Description             | Required    | Default value
        Tnbr         | name of the field       | no          | T
        kappa        | optional name of kappa  | no          | unknown
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            regionCoupledTemperature;
        Tnbr            T;
        kappa           kappa;
        value           uniform 300;
    }
    \endverbatim

    Needs to be on underlying mapped(Wall)FvPatch.

    The patch thermal conductivity \c kappa is obtained from the region
    thermophysicalTransportModel so that this boundary condition can be applied
    directly to either fluid or solid regions. If no 
    thermophysicalTransportModel is found, kappa is found by name.

SourceFiles
    regionCoupledTemperatureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef regionCoupledTemperatureFvPatchScalarField_H
#define regionCoupledTemperatureFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "scalarField.H"
#ifdef OPENFOAMFOUNDATION
#include "thermophysicalTransportModel.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class regionCoupledTemperatureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class regionCoupledTemperatureFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private Data

        //- Name of field on the neighbour region
        const word TnbrName_;

        //- Name of effective thermal conductivity field
        const word kappaName_;
        
        //- Flag to enable of disable the coupling
        bool enabled_;
        
    // Protected Functions
    
        //- Return the effective conductivity
        tmp<scalarField> kappa() const;
        
public:

    //- Runtime type information
    TypeName("regionCoupledTemperature");


    // Constructors

        //- Construct from patch and internal field
        regionCoupledTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );
        
        //- Construct from patch, internal field and dictionary
        regionCoupledTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  regionCoupledTemperatureFvPatchScalarField onto a
        //  new patch
        regionCoupledTemperatureFvPatchScalarField
        (
            const regionCoupledTemperatureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        regionCoupledTemperatureFvPatchScalarField
        (
            const regionCoupledTemperatureFvPatchScalarField&,
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
                new regionCoupledTemperatureFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }

    // Member Functions

        //- Return the thermophysical transport model
        const thermophysicalTransportModel& thermo() const;

        //- Is the coupling enabled
        inline bool enabled() const
        {
            return enabled_;
        }
        
        //- Enabled the coupled solution on this patchField
        inline void enable()
        {
            enabled_ = true;
        }
        
        //- Disable the coupled solution on this patchField
        inline void disable()
        {
            enabled_ = false;
        }
        
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
