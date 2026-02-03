/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2023 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "discontinuousCyclicAMIPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{


defineTypeNameAndDebug(discontinuousCyclicAMIPolyPatch, 0);

addToRunTimeSelectionTable(polyPatch, discontinuousCyclicAMIPolyPatch, word);
addToRunTimeSelectionTable(polyPatch, discontinuousCyclicAMIPolyPatch, dictionary);



// Constructors

discontinuousCyclicAMIPolyPatch::discontinuousCyclicAMIPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType,
    const transformType transform,
    const word& defaultAMIMethod
)
:
    cyclicAMIPolyPatch(name, size, start, index, bm, patchType, transform, defaultAMIMethod),
    duplicate_(false)
{}


discontinuousCyclicAMIPolyPatch::discontinuousCyclicAMIPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType,
    const word& defaultAMIMethod
)
:
    cyclicAMIPolyPatch(name, dict, index, bm, patchType, defaultAMIMethod),
    duplicate_(false)
{
    duplicate_ = dict.getOrDefault<bool>("duplicate", false);
}


discontinuousCyclicAMIPolyPatch::discontinuousCyclicAMIPolyPatch
(
    const discontinuousCyclicAMIPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    cyclicAMIPolyPatch(pp, bm),
    duplicate_(false)
{}


discontinuousCyclicAMIPolyPatch::discontinuousCyclicAMIPolyPatch
(
    const discontinuousCyclicAMIPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart,
    const word& nbrPatchName
)
:
    cyclicAMIPolyPatch(pp, bm, index, newSize, newStart, nbrPatchName),
    duplicate_(false)
{}


discontinuousCyclicAMIPolyPatch::discontinuousCyclicAMIPolyPatch
(
    const discontinuousCyclicAMIPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    cyclicAMIPolyPatch(pp, bm, index, mapAddressing, newStart),
    duplicate_(false)
{}


// Clone functions

autoPtr<polyPatch> discontinuousCyclicAMIPolyPatch::clone(const polyBoundaryMesh& bm) const
{
    return autoPtr<polyPatch>
    (
        new discontinuousCyclicAMIPolyPatch(*this, bm)
    );
}


autoPtr<polyPatch> discontinuousCyclicAMIPolyPatch::clone
(
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
) const
{
    return autoPtr<polyPatch>
    (
        new discontinuousCyclicAMIPolyPatch(*this, bm, index, newSize, newStart, neighbPatchName())
    );
}


autoPtr<polyPatch> discontinuousCyclicAMIPolyPatch::clone
(
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
) const
{
    return autoPtr<polyPatch>
    (
        new discontinuousCyclicAMIPolyPatch(*this, bm, index, mapAddressing, newStart)
    );
}

} // End namespace Foam


// ************************************************************************* //
