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
    Foam::blackBodyRadiationFvPatchScalarField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides a simplified black body radiation condition
    calculated as:

    \f[
        T_{mat} = \sqrt[4]{\frac{q''}{\sigma \epsilon} + T_a^4}
    \f]

    The heat flux is approximated using the Fourier law with the cells at the
    boundary:

    \f[
        q'' = \kappa \frac{T_{cell} - T_{patch face}}{d}
    \f]

    where
    \vartable
        T_{mat}     | Temperature of the material [K]
        q''         | Heat flux field [W/m2]
        \sigma      | Stefan-Boltzmann constant (5.670373e-8 W/(m2 K4))
        \epsilon    | Emissivity of the material surface [-]
        T_a         | Ambient temperature [K]
        T_{cell}    | Temperature of the cell [K]
        T_{surface} | Temperature of the patch face [K]
        \kappa      | Thermal conductivity of the material [W/(m K)]
        d           | Distance from the cell to patch face [m]
    \endtable

Usage
    \table
        Property   | Description                          | Req'd | Default
        emissivity | Emissivity of the material surface   | yes |
        kappa      | Thermal conductivity of the material | yes |
        Ta         | Ambient temperature                  | yes |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            blackBodyRadiation;
        emissivity      0.78;
        kappa           10;
        Ta              10;
        value           uniform 300; // Surface temperature
    }
    \endverbatim

    If used in another solver (e.g laplacianFoam), add the following line in
    system/controlDict:
    \verbatim
    libs ("libFFSSrc2.so");
    \endverbatim

    Author of this file (and associated .C or included .H files):
    Thomas Guilbaud <thomas.guilbaud@epfl.ch>, EPFL (Switzerland)

SourceFiles
    blackBodyRadiationFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef blackBodyRadiationFvPatchScalarField_H
#define blackBodyRadiationFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
     Class blackBodyRadiationFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class blackBodyRadiationFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private Data

        //- Emissivity of the medium
        scalar emissivity_;

        //- Thermal conductivity of the medium
        scalar kappa_;

        //- Ambient temperature
        scalar Ta_;


public:

    //- Runtime type information
    TypeName("blackBodyRadiation");


    // Constructors

        //- Construct from patch and internal field
        blackBodyRadiationFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        blackBodyRadiationFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        blackBodyRadiationFvPatchScalarField
        (
            const blackBodyRadiationFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy construct
        blackBodyRadiationFvPatchScalarField
        (
            const blackBodyRadiationFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new blackBodyRadiationFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        blackBodyRadiationFvPatchScalarField
        (
            const blackBodyRadiationFvPatchScalarField&,
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
                new blackBodyRadiationFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
