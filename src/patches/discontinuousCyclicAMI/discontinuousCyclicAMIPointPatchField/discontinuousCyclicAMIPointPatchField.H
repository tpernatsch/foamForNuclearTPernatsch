/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::cyclicAMIPointPatchField

Description
    Cyclic AMI front and back plane patch field

SourceFiles
    cyclicAMIPointPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_discontinuousCyclicAMIPointPatchField_H
#define Foam_discontinuousCyclicAMIPointPatchField_H

#include "cyclicAMIPointPatchField.H"
#include "discontinuousCyclicAMIPointPatch.H"
#include "primitivePatchInterpolation.H"

namespace Foam
{

template<class Type>
class discontinuousCyclicAMIPointPatchField
:
    public cyclicAMIPointPatchField<Type>
{
public:

    //- Runtime type information
    TypeName(discontinuousCyclicAMIPointPatch::typeName_());

    // Constructors

        discontinuousCyclicAMIPointPatchField
        (
            const pointPatch& patch,
            const DimensionedField<Type, pointMesh>& iF
        );

        discontinuousCyclicAMIPointPatchField
        (
            const pointPatch& patch,
            const DimensionedField<Type, pointMesh>& iF,
            const dictionary& dict
        );

        discontinuousCyclicAMIPointPatchField
        (
            const discontinuousCyclicAMIPointPatchField<Type>& ptf,
            const pointPatch& patch,
            const DimensionedField<Type, pointMesh>& iF,
            const pointPatchFieldMapper& mapper
        );

        discontinuousCyclicAMIPointPatchField
        (
            const discontinuousCyclicAMIPointPatchField<Type>& ptf,
            const DimensionedField<Type, pointMesh>& iF
        );

    // Member Functions

        //- Return the constraint type (disabled for discontinuous)
        virtual const word& constraintType() const
        {
            return word::null;
        }
};

} // End namespace Foam


#ifdef NoRepository
    #include "discontinuousCyclicAMIPointPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //