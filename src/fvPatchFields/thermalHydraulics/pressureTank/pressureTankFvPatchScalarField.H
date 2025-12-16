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
    Foam::pressureTankFvPatchScalarField

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
    pressureTankFvPatchScalarField.C
\*---------------------------------------------------------------------------*/


#ifndef pressureTankFvPatchScalarField_H
#define pressureTankFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class pressureTankFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/


class pressureTankFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{

protected:

    // Protected data

    scalar      level_;
    scalar      area_;
    scalar      p0_;
    scalar      currTimeIndex_;
    scalar      totalMass_;
    scalar      deltaM_;
    

public:

    //- Runtime type information
    TypeName("pressureTank");


    // Constructors
 
        //- Construct from patch and internal field
        pressureTankFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );
 
        //- Construct from patch, internal field and dictionary
        pressureTankFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );
 
        //- Construct by mapping given pressureTankFvPatchScalarField
        //  onto a new patch
        pressureTankFvPatchScalarField
        (
            const pressureTankFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );
 
        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new pressureTankFvPatchScalarField(*this)
            );
        }
 
        //- Construct as copy
        pressureTankFvPatchScalarField
        (
            const pressureTankFvPatchScalarField&
        );
 
        //- Construct as copy setting internal field reference
        pressureTankFvPatchScalarField
        (
            const pressureTankFvPatchScalarField&,
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
                new pressureTankFvPatchScalarField(*this, iF)
            );
        }

    // Member functions

        // Access

           virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #ifdef NoRepository
//     #include "pressureTankFvPatchScalarField.C"
// #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
