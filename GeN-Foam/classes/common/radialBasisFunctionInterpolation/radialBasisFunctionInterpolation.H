/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2306                                                  |
|    Copyright 2011-2016 OpenFOAM Foundation, 2017-2022 OpenCFD Ltd.         |
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
    Foam::radialBasisFunctionInterpolation

Description
    Base file for list conversion.

\mainauthor 
    Thomas Guilbaud (thomas.guilbaud@epfl.ch)

SourceFiles
    radialBasisFunctionInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef radialBasisFunctionInterpolation_H
#define radialBasisFunctionInterpolation_H

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radialBasisFunctionInterpolation
{

//- Weight calculation for the Radial Basis Function interpolation
//  https://en.wikipedia.org/wiki/Radial_basis_function_interpolation
scalarList solveGaussianRadialBasisFunction
(
    const scalarList xList,
    const scalarList yList,
    const scalarList zList,
    const scalarList vList,
    const scalarList eps,
    SquareMatrix<scalar>& invRBFmatrix
);

//- Compute the interpolation value based on the Radial Basis Function weights
//  https://en.wikipedia.org/wiki/Radial_basis_function_interpolation
scalar gaussianRadialBasisFunction
(
    const scalarList w,
    const scalarList xList,
    const scalarList yList,
    const scalarList zList,
    const scalar x,
    const scalar y,
    const scalar z,
    const scalarList eps
);

//- Polyharmonic spline RBF solver
scalarList solvePolyharmonicSpline
(
    const scalarList xList,
    const scalarList yList,
    const scalarList zList,
    const scalarList vList,
    SquareMatrix<scalar>& invRBFmatrix
);

//- Polyharmonic spline RBF
scalar polyharmonicSpline
(
    const scalarList w,
    const scalarList xList,
    const scalarList yList,
    const scalarList zList,
    const scalar x,
    const scalar y,
    const scalar z
);

//- Experimental RBF solver using derivative information (XY plane only)
scalarList solvePolyharmonicSplineDerivative
(
    const scalarList xList,
    const scalarList yList,
    const scalarList zList,
    const scalarList vList,
    const scalarList dvxList,
    const scalarList dvyList
);

//- Experimental RBF using derivative information (XY plane only)
scalar polyharmonicSplineDerivative
(
    const scalarList w,
    const scalarList xList,
    const scalarList yList,
    const scalarList zList,
    const scalar x,
    const scalar y,
    const scalar z
);

//- Experimental RBF solver imposing an integral under the interpolation field
scalarList solvePolyharmonicSplineIntegral
(
    const scalarList xList,
    const scalarList yList,
    const scalarList zList,
    const scalarList vList,
    const scalar totalIntegral,
    const labelList& regionCells,
    const fvMesh& mesh
);

//- Experimental RBF imposing an integral under the interpolation field
scalar polyharmonicSplineIntegral
(
    const scalarList w,
    const scalarList xList,
    const scalarList yList,
    const scalarList zList,
    const scalar x,
    const scalar y,
    const scalar z
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radialBasisFunctionInterpolation

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
