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
    Foam::convectiveFromFluidFvPatchScalarField

Description

    This class implements a mixed boundary condition which is meant
    to be used in combination with a fluid model for multi-scale 
    simulations. The patch creates two additional fields:
    - htc: the heat transfer coefficient
    - Tfluid: the fluid temperature

    These additonal fields are created so that the user can   
    map the fluid paramters onto the solid patch using the
    volumetric mapping capabilties of foamForNuclear. Once
    the fields are mapped, they are used to set a mixed
    boundary condition on the solid patch that conserves
    the heat flux between the solid and the fluid. 

    If adjancent solid-fluid simulations are intended,
    then other patches or the foamForNuclear CHTLoop
    features should be used.

    Author of this file (and associated .C or included .H files):
    Giovanni Nervi <giovanni.nervi@epfl.ch>, EPFL (Switzerland)

SourceFiles
    convectiveFromFluidFvPatchScalarField.H

\*---------------------------------------------------------------------------*/


#ifndef convectiveFromFluidFvPatchScalarField_H
#define convectiveFromFluidFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "volFields.H"
#include "surfaceFields.H"

namespace Foam
{

class convectiveFromFluidFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data

    // conductivity field name (solid property)
    word kappaName_;

    // Htc and fluid T field names (created/owned here)
    word htcName_;
    word TfluidName_;

    volScalarField TFluid_;
    volScalarField htc_;

public:

    //- Runtime type information
    TypeName("convectiveFromFluid");


    // Constructors

    //- Construct from patch and internal field
    convectiveFromFluidFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    //- Construct from patch, internal field and dictionary
    convectiveFromFluidFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    //- Construct as copy with new internal field
    convectiveFromFluidFvPatchScalarField
    (
        const convectiveFromFluidFvPatchScalarField&,
        const DimensionedField<scalar, volMesh>&
    );


    //- Construct by mapping
    convectiveFromFluidFvPatchScalarField
    (
        const convectiveFromFluidFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    //- Clone
    virtual tmp<fvPatchScalarField> clone() const
    {
        return tmp<fvPatchScalarField>
        (
            new convectiveFromFluidFvPatchScalarField(*this)
        );
    }

    //- Clone with new internal field
    virtual tmp<fvPatchScalarField> clone
    (
        const DimensionedField<scalar, volMesh>& iF
    ) const
    {
        return tmp<fvPatchScalarField>
        (
            new convectiveFromFluidFvPatchScalarField(*this, iF)
        );
    }


    // Member functions

    //- Update coefficients
    virtual void updateCoeffs();

    //- Write
    virtual void write(Ostream&) const;
};

} // End namespace Foam

#endif
