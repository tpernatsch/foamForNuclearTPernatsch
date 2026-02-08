/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "regionCoupledOFFBEATPolyPatch.H"
#include "polyMesh.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(regionCoupledOFFBEATPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, regionCoupledOFFBEATPolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, regionCoupledOFFBEATPolyPatch, dictionary);
}


// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * //

Foam::regionCoupledOFFBEATPolyPatch::regionCoupledOFFBEATPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    polyPatch(name, size, start, index, bm, patchType),
    regionCoupledBaseOFFBEAT(static_cast<const polyPatch&>(*this))
{}


Foam::regionCoupledOFFBEATPolyPatch::regionCoupledOFFBEATPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    polyPatch(name, dict, index, bm, patchType),
    regionCoupledBaseOFFBEAT(*this, dict)
{}


Foam::regionCoupledOFFBEATPolyPatch::regionCoupledOFFBEATPolyPatch
(
    const regionCoupledOFFBEATPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    polyPatch(pp, bm),
    regionCoupledBaseOFFBEAT(*this, pp)
{}


Foam::regionCoupledOFFBEATPolyPatch::regionCoupledOFFBEATPolyPatch
(
    const regionCoupledOFFBEATPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    polyPatch(pp, bm, index, newSize, newStart),
    regionCoupledBaseOFFBEAT(*this, pp)
{}


Foam::regionCoupledOFFBEATPolyPatch::regionCoupledOFFBEATPolyPatch
(
    const regionCoupledOFFBEATPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    polyPatch(pp, bm, index, mapAddressing, newStart),
    regionCoupledBaseOFFBEAT(*this, pp)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionCoupledOFFBEATPolyPatch::~regionCoupledOFFBEATPolyPatch()
{
    regionCoupledBaseOFFBEAT::clearGeom();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionCoupledOFFBEATPolyPatch::initGeometry(PstreamBuffers& pBufs)
{
#ifdef OPENFOAMFOUNDATION    
    polyPatch::initCalcGeometry(pBufs);
#elif OPENFOAMESI    
    polyPatch::initGeometry(pBufs);
#endif    
}


void Foam::regionCoupledOFFBEATPolyPatch::initMovePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    polyPatch::initMovePoints(pBufs, p);
}


void Foam::regionCoupledOFFBEATPolyPatch::movePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    polyPatch::movePoints(pBufs, p);
    regionCoupledBaseOFFBEAT::clearGeom();
}


void Foam::regionCoupledOFFBEATPolyPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    polyPatch::initUpdateMesh(pBufs);
}


void Foam::regionCoupledOFFBEATPolyPatch::updateMesh(PstreamBuffers& pBufs)
{
    polyPatch::updateMesh(pBufs);
    regionCoupledBaseOFFBEAT::clearGeom();
}


void Foam::regionCoupledOFFBEATPolyPatch::write(Ostream& os) const
{
    polyPatch::write(os);
    regionCoupledBaseOFFBEAT::write(os);
}


// ************************************************************************* //
