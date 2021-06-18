/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "regimeDomain2D.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regimeMapModels
{
    defineTypeNameAndDebug(regimeDomain2D, 0);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regimeMapModels::regimeDomain2D::regimeDomain2D()
:
    id_(-1)
{}


Foam::regimeMapModels::regimeDomain2D::regimeDomain2D
(
    label id,
    List<const Vector2D<scalar>*> pointPtrs
)
:
    id_(id),
    pointPtrs_(pointPtrs),
    boundaries_(0),
    internalBoundaryPtrs_(0)
{
    //Info << "    constructing regimeDomain2D ID " << id_ << endl;
    int N(pointPtrs_.size());
    forAll(pointPtrs_, i)
    {
        const Vector2D<scalar>* p0Ptr(nullptr);
        const Vector2D<scalar>* p1Ptr(nullptr);
        if (i < N-1)
        {
            p0Ptr = pointPtrs_[i];
            p1Ptr = pointPtrs_[i+1];
        }
        else
        {
            p0Ptr = pointPtrs_[N-1];
            p1Ptr = pointPtrs_[0];
        }
        boundaries_.append
        (
            regimeBoundary2D
            (
                p0Ptr,
                p1Ptr
            )
        );
        //Info<< "        added boundary " << boundaries_[i].p0() << " " 
        //    << boundaries_[i].p1() << endl;
    }

    //- I should implement sanity-checks to ensure that the user-provided 
    //  points from a convex polygon, with a counter-clockwise point ordering.
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::Tuple2<Foam::Vector2D<Foam::scalar>, Foam::Vector2D<Foam::scalar>> 
Foam::regimeMapModels::regimeDomain2D::boundingBox() const
{
    const Vector2D<scalar>& p0(*pointPtrs_[0]);
    scalar minX=p0[0], minY=p0[1], maxX=p0[0], maxY=p0[1];

    forAll(pointPtrs_, i)
    {
        const Vector2D<scalar>& p(*pointPtrs_[i]);
        minX = min(minX, p[0]);
        maxX = max(maxX, p[0]);
        minY = min(minY, p[1]);
        maxY = max(maxY, p[1]);
    }

    return 
        Tuple2<Vector2D<scalar>, Vector2D<scalar>>
        (
            Vector2D<scalar>(minX, minY),
            Vector2D<scalar>(maxX, maxY)
        );
}


void Foam::regimeMapModels::regimeDomain2D::setInternalBoundaries()
{
    forAll(boundaries_, i)
    {
        const regimeBoundary2D& b(boundaries_[i]);
        if (b.nbrId() != -1)
            internalBoundaryPtrs_.append(&b);
    }
}


bool Foam::regimeMapModels::regimeDomain2D::containsPoint
(
    const Vector2D<scalar>& p
) const
{
    //- The logic here is to perform a cross-product (wherein the thrid cmpt of
    //  the involved vectors is always 0, as these objects are 2-D) between
    //  each of the vectors that can be constructed between each boundary
    //  starting poing p0 and the point under exam p with the corresponding
    //  p0-p1 vector for each boundary, for all boundaries. If the sign of 
    //  these products is consistently positive for all boundaries, then the
    //  point lies in the domain. Note that this ONLY works if BOTH these
    //  conditions are satisfied:
    //  
    //  -   the point list used to construct the regimeDomain is provided in a
    //      counter-clockwise (i.e. right-hand-rule) order;
    //
    //  -   the regime domain is convex.
    //
    //  I suppose it would not be too difficult to generalize this to points
    //  provided in a clockwise order, and for non-covex regime domains, yet
    //  for now I'll leave it at this, it is already quite flexible
    forAll(boundaries_, i)
    {
        const regimeBoundary2D& b(boundaries_[i]);
        const Vector2D<scalar> v01(b.v());
        Vector2D<scalar> v0p(p-b.p0());
        scalar x(v0p[1]*v01[0] - v0p[0]*v01[1]);
        if (x < 0)
            return false;
    }
    return true;
}


// ************************************************************************* //
