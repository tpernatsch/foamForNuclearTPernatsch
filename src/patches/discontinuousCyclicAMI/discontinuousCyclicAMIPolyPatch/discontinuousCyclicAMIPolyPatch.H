/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2023 OpenCFD Ltd.
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
    Foam::discontinuousdiscontinuousCyclicAMIPolyPatch

Description
    Cyclic patch for Arbitrary Mesh Interface (AMI)

    Includes provision for updating the patch topology to enforce a 1-to-1
    face match across the interface, based on the \c createAMIFaces flag.

    The manipulations are based on the reference:

    \verbatim
    H.J. Aguerre, S. Márquez Damián, J.M. Gimenez, N.M.Nigro, Conservative
    handling of arbitrary non-conformal interfaces using an efficient
    supermesh, Journal of Computational Physics 335(15) 21-49. 2017.
    https://doi.org/10.1016/j.jcp.2017.01.018.
    \endverbatim

SourceFiles
    discontinuousCyclicAMIPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_discontinuousCyclicAMIPolyPatch_H
#define Foam_discontinuousCyclicAMIPolyPatch_H

#include "cyclicAMIPolyPatch.H"

namespace Foam
{

class discontinuousCyclicAMIPolyPatch
:
    public cyclicAMIPolyPatch
{

    protected:
    
        bool duplicate_;
public:

    //- Runtime type information
    TypeName("discontinuousCyclicAMI");

    // Constructors

    //- Construct from (base coupled patch) components
    discontinuousCyclicAMIPolyPatch
    (
        const word& name,
        const label size,
        const label start,
        const label index,
        const polyBoundaryMesh& bm,
        const word& patchType,
        const transformType transform = UNKNOWN,
        const word& defaultAMIMethod = faceAreaWeightAMI::typeName
    );

    //- Construct from dictionary
    discontinuousCyclicAMIPolyPatch
    (
        const word& name,
        const dictionary& dict,
        const label index,
        const polyBoundaryMesh& bm,
        const word& patchType,
        const word& defaultAMIMethod = faceAreaWeightAMI::typeName
    );

    //- Construct as copy, resetting the boundary mesh
    discontinuousCyclicAMIPolyPatch
    (
        const discontinuousCyclicAMIPolyPatch& pp,
        const polyBoundaryMesh& bm
    );

    //- Construct given the original patch and resetting the face list and bm
    discontinuousCyclicAMIPolyPatch
    (
        const discontinuousCyclicAMIPolyPatch& pp,
        const polyBoundaryMesh& bm,
        const label index,
        const label newSize,
        const label newStart,
        const word& nbrPatchName
    );

    //- Construct given the original patch and a map
    discontinuousCyclicAMIPolyPatch
    (
        const discontinuousCyclicAMIPolyPatch& pp,
        const polyBoundaryMesh& bm,
        const label index,
        const labelUList& mapAddressing,
        const label newStart
    );

    // Cloning (preserve derived type)

    virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const;

    virtual autoPtr<polyPatch> clone
    (
        const polyBoundaryMesh& bm,
        const label index,
        const label newSize,
        const label newStart
    ) const;

    virtual autoPtr<polyPatch> clone
    (
        const polyBoundaryMesh& bm,
        const label index,
        const labelUList& mapAddressing,
        const label newStart
    ) const;

    virtual bool duplicate() const
    {
        return duplicate_;
    }
};

} // End namespace Foam


// ************************************************************************* //

#endif

// ************************************************************************* //
