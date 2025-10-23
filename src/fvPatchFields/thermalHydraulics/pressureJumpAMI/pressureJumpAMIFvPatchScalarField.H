/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2506                                                  |
|    Copyright 2011-2016 OpenFOAM Foundation, 2017-2025 OpenCFD Ltd.          |
-------------------------------------------------------------------------------
License
    This file is part of GeN-Foam.

    GeN-Foam is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    GeN-Foam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    This offering is not approved or endorsed by the OpenFOAM Foundation nor
    OpenCFD Limited, producer and distributor of the OpenFOAM(R)software via
    www.openfoam.com, and owner of the OPENFOAM(R) and OpenCFD(R) trademarks.

    This particular snippet of code is developed according to the developer's
    knowledge and experience in OpenFOAM. The users should be aware that
    there is a chance of bugs in the code, though we've thoroughly test it.
    The source code may not be in the OpenFOAM coding style, and it might not
    be making use of inheritance of classes to full extent.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::pressureJumpAMIFvPatchScalarField

Group
    grpCoupledBoundaryConditions

Author
    Giovanni Nervi (giovanni.nervi@epfl.ch)

Description
    This boundary condition imposes a pressure jump across non-conformal
    cyclic path pairs using an arbitraryMeshInterface (AMI). The total jump
    is composed of four contributions:

    1) A user-defined static pressure source (in Pa);
    2) A Bernoulli term accounting for velocity-induced pressure changes;
    3) A pressure loss term using a user-defined loss coefficient;
    4) A pump contribution defined by a characteristic curve with
       user-specified coefficients.

Usage
    \table
        Property         | Description                                                  | Required | Default
        patchType        | underlying patch type (should be 'cyclicAMI')                | yes      |
        pSource          | static pressure source in Pa                                 | no       | 0
        bernoulli        | switch to enable Bernoulli contribution                      | no       | true
        lossCoeff        | loss coefficient for pressure drop                           | no       | 0
        pumpCoeffs       | list of coefficients defining the pump characteristic curve  | no       |(0)
        pumpSpeedTable   | Table containing pump rotational velocity vs time            | no       |
        nominalVelocity  | Speed at which nominal characteristic curve is defined       | no       |
    \endtable

    Example:
    \verbatim
    <patchName>
    {
        type            bernoulliJumpAMI;
        patchType       discontinuousCyclicAMI;

        // 1) Static pressure contribution in Pa
        pSource         uniform 5000;

        // 2) Bernoulli contribution
        bernoulli       on;

        // 3) Loss coefficient
        lossCoeff       uniform 0.2;

        // 4) Pump characteristic curve coefficients (example: a0 + a1*Q + a2*Q^2)
        pumpCoeffs      (10000 -500 20);
    }
    \endverbatim

Note
    The underlying patchType must be 'cyclicAMI'. The total applied jump is
    the sum of the user-defined pressure source, the Bernoulli contribution,
    the loss-induced drop, and the pump characteristic term evaluated
    from the supplied coefficients.

See also
    Foam::jumpDiscontinuousCyclicAMIFvPatchField

SourceFiles
    pressureJumpAMIFvPatchScalarField.C
\*---------------------------------------------------------------------------*/


#ifndef pressureJumpAMIFvPatchScalarField_H
#define pressureJumpAMIFvPatchScalarField_H

#include "jumpDiscontinuousCyclicAMIFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class pressureJumpAMIFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/


class pressureJumpAMIFvPatchScalarField
:
    public jumpDiscontinuousCyclicAMIFvPatchField<scalar>
{

protected:

    // Protected data

    scalarField                 jump_; // total jump applied: jump = bernoulli + source - losses + pump
    bool                        bernoulliCorrection_; // option to include bernoulli pressure correction
    scalarField                 source_; // user-input total source
    scalar                      pressureLossCoeff_; // K in delta_P = K * 0.5 * rho * magU^2
    scalarList                  pumpCoeffs_;  // pump characteristics coefficients 
    autoPtr<Function1<scalar>>  pumpVelocityTable_; // pump speed table to correct characteristics
    scalar                      nominalVelocity_;
    

public:

    //- Runtime type information
    TypeName("pressureJumpAMI");


    // Constructors

        //- Construct from patch and internal field
        pressureJumpAMIFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        pressureJumpAMIFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given pressureJumpAMIFvPatchScalarField onto a
        //  new patch
        pressureJumpAMIFvPatchScalarField
        (
            const pressureJumpAMIFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        pressureJumpAMIFvPatchScalarField
        (
            const pressureJumpAMIFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        pressureJumpAMIFvPatchScalarField
        (
            const pressureJumpAMIFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return fvPatchField<scalar>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return fvPatchField<scalar>::Clone(*this, iF);
        }


    // Member functions

        // Access

            //- Return the "jump" across the patch
            virtual tmp<Field<scalar>> jump() const;


        // Mapping functions

            virtual void updateCoeffs();

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap(const fvPatchField<scalar>&, const labelList&);


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #ifdef NoRepository
//     #include "pressureJumpAMIFvPatchScalarField.C"
// #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
