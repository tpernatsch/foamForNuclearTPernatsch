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
    Foam::discontinuousCyclicAMIFvsPatchField

Author
    Giovanni Nervi (giovanni.nervi@epfl.ch)

Description
    This boundary condition is a modification of the OpenFOAM cyclicAMIFvPatchField.
    This boundary conditions switches between a zeroGradient BC and a jump cyclicBC,
    depending on whether the patchField is being created from dictionary (i.e., read
    from file) or being created during the simulation. The two modes of the patchField
    allow to have discontinuity in certain fields (e.g., porosity, hydraulic diameter)
    which are not physically continuous across AMI interface, while still allowing
    other fields (e.g., velocity, pressure) to be coupled.

    The "value" entry is MUST_READ.

SourceFiles
    discontinuousCyclicAMIFvsPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_discontinuousCyclicAMIFvsPatchField_H
#define Foam_discontinuousCyclicAMIFvsPatchField_H

#include "discontinuousCyclicAMIFvPatch.H"
#include "coupledFvsPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class discontinuousCyclicAMIFvsPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class discontinuousCyclicAMIFvsPatchField
:
    public coupledFvsPatchField<Type>
{
    // Private Data

        //- Local reference cast into the cyclic patch
        const discontinuousCyclicAMIFvPatch& discontinuousCyclicAMIPatch_;


public:

    //- Runtime type information
    TypeName(discontinuousCyclicAMIFvPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        discontinuousCyclicAMIFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Construct from patch, internal field and dictionary.
        //- The "value" entry is MUST_READ.
        discontinuousCyclicAMIFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const dictionary&
        );

        //- Construct by mapping given discontinuousCyclicAMIFvsPatchField onto a new patch
        discontinuousCyclicAMIFvsPatchField
        (
            const discontinuousCyclicAMIFvsPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy -- From OpenFOAM v2606, this needs to be deleted
        discontinuousCyclicAMIFvsPatchField
        (
            const discontinuousCyclicAMIFvsPatchField<Type>&
        ) = delete;

        //- Construct as copy setting internal field reference
        discontinuousCyclicAMIFvsPatchField
        (
            const discontinuousCyclicAMIFvsPatchField<Type>&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Return clone -- Changed to compile in v2606
        virtual tmp<fvsPatchField<Type>> clone() const
        {
            // Use the constructor that takes (ptf, internalField)
            return fvsPatchField<Type>::Clone(*this, this->internalField());
        }


        //- Clone with an internal field reference
        virtual tmp<fvsPatchField<Type>> clone
        (
            const DimensionedField<Type, surfaceMesh>& iF
        ) const
        {
            return fvsPatchField<Type>::Clone(*this, iF);
        }


    // Member Functions

        //- Return true if running parallel
        virtual bool coupled() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "discontinuousCyclicAMIFvsPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
