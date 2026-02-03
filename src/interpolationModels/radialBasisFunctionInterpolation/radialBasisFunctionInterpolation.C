/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2512                                                  |
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

\*---------------------------------------------------------------------------*/

#include "radialBasisFunctionInterpolation.H"
#include "LUscalarMatrix.H"


namespace Foam
{

namespace radialBasisFunctionInterpolation
{


scalarList solveGaussianRadialBasisFunction
(
    const scalarList xList,
    const scalarList yList,
    const scalarList zList,
    const scalarList vList,
    const scalarList eps,
    SquareMatrix<scalar>& invRBFmatrix
)
{
    const label nx(xList.size());

    scalarList w(nx, 0.0);

    if (invRBFmatrix.m() != nx)
    {
        invRBFmatrix.resize(nx);

        SquareMatrix<scalar> A(nx, 0.0);

        forAll(xList, i)
        {
            forAll(xList, j)
            {
                A[i][j] = exp(
                    -eps[0] * sqr(xList[i]-xList[j])
                    -eps[1] * sqr(yList[i]-yList[j])
                    -eps[2] * sqr(zList[i]-zList[j])
                );
            }
        }

        LUscalarMatrix Atemp(A);

        Atemp.inv(invRBFmatrix);
    }

    w = invRBFmatrix * vList;

    return(w);
}


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
)
{
    scalar res(0);
    forAll(w, i)
    {
        res += w[i] * exp(
            -eps[0] * sqr(x-xList[i])
            -eps[1] * sqr(y-yList[i])
            -eps[2] * sqr(z-zList[i])
        );
    }
    return(res);
}


scalar polyharmonicSplineFunction
(
    const scalar rSquare,
    const label mode
)
{
    switch (mode) {
        case 1:
            return(sqrt(rSquare));
        case 2:
            return(rSquare * 0.5 * log(rSquare));
        case 3:
            return(sqrt(rSquare) * rSquare);
        case 4:
            return(sqr(rSquare) * 0.5 * log(rSquare));
        default:
            Info<< "Polyharmonic spline mode " << mode << " not in range [1; 4], return 0"
                << endl;
            break;
    }
    return(0);
}


scalarList solvePolyharmonicSpline
(
    const scalarList xList,
    const scalarList yList,
    const scalarList zList,
    const scalarList vList,
    SquareMatrix<scalar>& invRBFmatrix,
    const label mode
)
{
    const label nx(xList.size());

    scalarList w(nx+4, 0.0);

    if (invRBFmatrix.m() != nx+4)
    {
        invRBFmatrix.resize(nx+4);

        SquareMatrix<scalar> A(nx+4, 0.0);

        scalar rSquare(0);
        forAll(xList, i)
        {
            forAll(xList, j)
            {
                if (i != j)
                {
                    rSquare = sqr(xList[i]-xList[j])
                            + sqr(yList[i]-yList[j])
                            + sqr(zList[i]-zList[j]);
                    A[i][j] = polyharmonicSplineFunction(rSquare, mode);
                }
                else
                {
                    A[i][j] = 0.0;
                }
            }
            // Polynomial correction
            A[i][nx] = 1.0;
            A[i][nx+1] = xList[i];
            A[i][nx+2] = yList[i];
            A[i][nx+3] = zList[i];
            A[nx][i] = 1.0;
            A[nx+1][i] = xList[i];
            A[nx+2][i] = yList[i];
            A[nx+3][i] = zList[i];
        }

        // Inverse the matrix once and store it for later iterations
        LUscalarMatrix Atemp(A);

        Atemp.inv(invRBFmatrix);
    }

    scalarList vListTemp(w.size(), 0.0);
    forAll(vList, paramI)
    {
        vListTemp[paramI] = vList[paramI];
    }

    w = invRBFmatrix * vListTemp;

    return(w);
}

scalarList solvePolyharmonicSpline
(
    const List<scalarList>& xList,
    const scalarList& vList,
    SquareMatrix<scalar>& invRBFmatrix,
    const label mode
)
{
    const label nx(xList.first().size());
    const label nParameters(xList.size());

    scalarList w(nx+nParameters+1, 0.0);

    if (invRBFmatrix.m() != w.size())
    {
        invRBFmatrix.resize(w.size());

        SquareMatrix<scalar> A(w.size(), 0.0);

        scalar rSquare(0);
        forAll(xList.first(), dataI)
        {
            forAll(xList.first(), dataJ)
            {
                if (dataI != dataJ)
                {
                    rSquare = 0.0;
                    forAll(xList, paramI)
                    {
                        rSquare += sqr(xList[paramI][dataI] - xList[paramI][dataJ]);
                    }
                    scalar poly(polyharmonicSplineFunction(rSquare, mode));
                    A[dataI][dataJ] = poly;
                    A[dataJ][dataI] = poly;
                }
                else
                {
                    A[dataI][dataJ] = 0.0;
                }
            }
            // Polynomial correction
            A[dataI][nx] = 1.0;
            A[nx][dataI] = 1.0;
            forAll(xList, paramI)
            {
                A[dataI][nx+paramI+1] = xList[paramI][dataI];
                A[nx+paramI+1][dataI] = xList[paramI][dataI];
            }
        }

        // Inverse the matrix once and store it for later iterations
        LUscalarMatrix Atemp(A);

        Atemp.inv(invRBFmatrix);
    }

    scalarList vListTemp(w.size(), 0.0);
    forAll(vList, paramI)
    {
        vListTemp[paramI] = vList[paramI];
    }

    w = invRBFmatrix * vListTemp;

    return(w);
}


scalar polyharmonicSpline
(
    const scalarList w,
    const scalarList xList,
    const scalarList yList,
    const scalarList zList,
    const scalar x,
    const scalar y,
    const scalar z,
    const label mode
)
{
    const label nx(xList.size());
    scalar res(0);
    scalar rSquare(0);
    forAll(xList, i)
    {
        rSquare = sqr(x - xList[i])
                + sqr(y - yList[i])
                + sqr(z - zList[i]);
        res += w[i] * polyharmonicSplineFunction(rSquare, mode);
    }
    res += w[nx] + w[nx+1]*x + w[nx+2]*y + w[nx+3]*z;
    return(res);
}

scalar polyharmonicSpline
(
    const scalarList& w,
    const List<scalarList>& xList,
    const scalarList& xInput,
    const label mode
)
{
    const label nx(xList.first().size());

    // Polynomial Correction
    scalar res(w[nx]);
    forAll(xInput, paramI)
    {
        res += w[nx+paramI+1] * xInput[paramI];
    }

    scalar rSquare(0);
    forAll(xList.first(), dataI)
    {
        rSquare = 0.0;
        forAll(xInput, paramI)
        {
            rSquare += sqr(xInput[paramI] - xList[paramI][dataI]);
        }
        if (rSquare > 0)
        {
            res += w[dataI] * polyharmonicSplineFunction(rSquare, mode);
        }
    }
    return(res);
}


scalarList solvePolyharmonicSplineDerivative
(
    const scalarList xList,
    const scalarList yList,
    const scalarList zList,
    const scalarList vList,
    const scalarList dvxList,
    const scalarList dvyList
)
{
    const label nx(xList.size());
    SquareMatrix<scalar> A(3*nx+4, 0.0);
    scalarList w(3*nx+4, 0.0);

    forAll(xList, i)
    {
        const scalar xi(xList[i]);
        const scalar yi(yList[i]);
        const scalar zi(zList[i]);
        forAll(xList, j)
        {
            const scalar xj(xList[j]);
            const scalar yj(yList[j]);
            const scalar dx(xi-xj);
            const scalar dy(yi-yj);
            const scalar r(sqrt(
                sqr(dx) + sqr(dy) + sqr(zi-zList[j])
            ));
            const scalar rd(sqrt(sqr(dx) + sqr(dy)));
            const scalar phir(sqr(r) * log(r));
            if (r > 0)
            {

                // wi with fi (A)
                A[i][j] = phir;
                // vix with fi (C)
                A[i][nx+j] = dx * phir;
                // viy with fi (H)
                A[i][2*nx+j] = dy * phir;
            }
            if (rd > 0)
            {
                const scalar dphix(2.0 * dx * log(rd) + dx);
                const scalar dphiy(2.0 * dy * log(rd) + dy);
                // wi with df/dx (D)
                A[nx+i][j] = dphix;
                // vix with df/dx (E)
                A[nx+i][nx+j] = phir + dx * dphix;
                // viy with df/dx (F)
                A[nx+i][2*nx+j] = dy * dphix;

                // wi with df/dy (I)
                A[2*nx+i][j] = dphiy;
                // vix with df/dy (K)
                A[2*nx+i][nx+j] = dx * dphiy;
                // viy with df/dy (J)
                A[2*nx+i][2*nx+j] = phir + dy * dphiy;
            }
        }
        // Polynomial correction (B)
        A[i][3*nx] = 1;
        A[i][3*nx+1] = xi;
        A[i][3*nx+2] = yi;
        A[i][3*nx+3] = zi;
        A[3*nx][i] = 1;
        A[3*nx+1][i] = xi;
        A[3*nx+2][i] = yi;
        A[3*nx+3][i] = zi;

        // Polynome for df/dx (G)
        A[nx+i][3*nx+1] = 1;
        // Polynome for df/dy (L)
        A[2*nx+i][3*nx+2] = 1;
    }

    scalarList vListTemp(vList);
    forAll(dvxList, i)
    {
        vListTemp.append(dvxList[i]);
    }
    forAll(dvyList, i)
    {
        vListTemp.append(dvyList[i]);
    }
    for (int i = 0; i < 4; i++)
    {
        vListTemp.append(0);
    }

    solve(w, A, vListTemp);

    return(w);
}


scalar polyharmonicSplineDerivative
(
    const scalarList w,
    const scalarList xList,
    const scalarList yList,
    const scalarList zList,
    const scalar x,
    const scalar y,
    const scalar z
)
{
    const label nx(xList.size());
    scalar res(0);
    scalar rSquare(0);
    forAll(xList, i)
    {
        rSquare = sqr(x-xList[i])
                + sqr(y-yList[i])
                + sqr(z-zList[i]);
        if (rSquare > 0)
        {
            res += (
                w[i] + w[i+nx] * (x-xList[i]) + w[i+2*nx] * (y-yList[i])
            ) * rSquare * 0.5 * log(rSquare);
        }
    }
    res += w[3*nx] + w[3*nx+1]*x + w[3*nx+2]*y + w[3*nx+3]*z;
    return(res > 0 ? res : 0);
}


scalarList solvePolyharmonicSplineIntegral
(
    const scalarList xList,
    const scalarList yList,
    const scalarList zList,
    const scalarList vList,
    const scalar totalIntegral,
    const labelList& regionCells,
    const fvMesh& mesh,
    SquareMatrix<scalar>& invRBFmatrix,
    const label mode
)
{
    const label nx(xList.size());

    scalarList w(nx+1, 0.0);

    if (invRBFmatrix.m() != nx+1)
    {
        invRBFmatrix.resize(nx+1);

        SquareMatrix<scalar> A(nx+1, 0.0);

        const scalarList& V(mesh.V());

        scalar rSquare(0);
        forAll(xList, i)
        {
            forAll(xList, j)
            {
                if (i != j)
                {
                    rSquare = sqr(xList[i] - xList[j])
                            + sqr(yList[i] - yList[j])
                            + sqr(zList[i] - zList[j]);
                    A[i][j] = polyharmonicSplineFunction(rSquare, mode);
                }
                else
                {
                    A[i][j] = 0.0;
                }
            }
            // Polynomial correction
            A[i][nx] = 1.0;

            scalar totalPhi(0), totalVolume(0);
            forAll(regionCells, celli)
            {
                rSquare = sqr(mesh.C().internalField()[celli].x() - xList[i])
                        + sqr(mesh.C().internalField()[celli].y() - yList[i])
                        + sqr(mesh.C().internalField()[celli].z() - zList[i]);

                totalPhi += polyharmonicSplineFunction(rSquare, mode) * V[celli];
                totalVolume += V[celli];
            }

            A[nx][i] = totalPhi;
            A[nx][nx] = totalVolume;
        }

        // Inverse the matrix once and store it for later iterations
        LUscalarMatrix Atemp(A);

        Atemp.inv(invRBFmatrix);
    }

    scalarList vListTemp(vList);
    vListTemp.append(totalIntegral);

    w = invRBFmatrix * vListTemp;

    return(w);
}


scalar polyharmonicSplineIntegral
(
    const scalarList w,
    const scalarList xList,
    const scalarList yList,
    const scalarList zList,
    const scalar x,
    const scalar y,
    const scalar z,
    const label mode
)
{
    const label nx(xList.size());
    scalar res(w[nx]);
    scalar rSquare(0);
    forAll(xList, i)
    {
        rSquare = sqr(x - xList[i])
                + sqr(y - yList[i])
                + sqr(z - zList[i]);
        res += w[i] * polyharmonicSplineFunction(rSquare, mode);
    }
    return(res);
}



void solveKriging
(
    const scalarList xList,
    const scalarList yList,
    const scalarList zList,
    const variogramType type,
    const scalar a,
    const scalar c,
    const scalar c0,
    SquareMatrix<scalar>& invRBFmatrix
)
{
    const label nx(xList.size());

    if (invRBFmatrix.m() != nx+1)
    {
        invRBFmatrix.resize(nx+1);

        SquareMatrix<scalar> A(nx+1, 0.0);

        scalar distance(0);
        switch (type) {
            case spherical:
                forAll(xList, i)
                {
                    forAll(xList, j)
                    {
                        distance = sqrt(
                            sqr(xList[i]-xList[j]) +
                            sqr(yList[i]-yList[j]) +
                            sqr(zList[i]-zList[j])
                        ) / a;
                        if (distance <= 1) {
                            A[i][j] = c0 + c*(1.5*distance - 0.5*pow(distance, 3.0));
                        } else {
                            A[i][j] = c0 + c;
                        }
                    }
                    A[i][nx] = 1.0;
                    A[nx][i] = 1.0;
                }
                break;
            case gaussian:
                forAll(xList, i)
                {
                    forAll(xList, j)
                    {
                        distance = sqr(xList[i]-xList[j]) +
                            sqr(yList[i]-yList[j]) +
                            sqr(zList[i]-zList[j]);
                        A[i][j] = c0 + c * (1.0 - exp(-distance / sqr(a)));
                    }
                    A[i][nx] = 1.0;
                    A[nx][i] = 1.0;
                }
                break;
            case exponential:
                forAll(xList, i)
                {
                    forAll(xList, j)
                    {
                        distance = sqrt(
                            sqr(xList[i]-xList[j]) +
                            sqr(yList[i]-yList[j]) +
                            sqr(zList[i]-zList[j])
                        );
                        A[i][j] = c0 + c * (1.0 - exp(-distance/a));
                    }
                    A[i][nx] = 1.0;
                    A[nx][i] = 1.0;
                }
                break;
        }

        LUscalarMatrix Atemp(A);

        Atemp.inv(invRBFmatrix);
    }
}


scalarList getVariogramVector
(
    const scalarList xList,
    const scalarList yList,
    const scalarList zList,
    const variogramType type,
    const scalar a,
    const scalar c,
    const scalar c0,
    const scalar x,
    const scalar y,
    const scalar z
)
{
    const label nx(xList.size());

    scalarList xVariogram(nx+1);
    scalar distance(0);
    switch (type) {
        case spherical:
            forAll(xList, i)
            {
                distance = sqrt(
                    sqr(xList[i]-x) + sqr(yList[i]-y) + sqr(zList[i]-z)
                ) / a;
                if (distance <= 1) {
                    xVariogram[i] = c0 + c*(1.5*distance - 0.5*pow(distance, 3.0));
                } else {
                    xVariogram[i] = c0 + c;
                }
            }
            break;
        case gaussian:
            forAll(xList, i)
            {
                distance = sqr(xList[i]-x) +  sqr(yList[i]-y) + sqr(zList[i]-z);
                xVariogram[i] = c0 + c * (1.0 - exp(-distance/sqr(a)));
            }
            break;
        case exponential:
            forAll(xList, i)
            {
                distance = sqrt(
                    sqr(xList[i]-x) +
                    sqr(yList[i]-y) +
                    sqr(zList[i]-z)
                );
                xVariogram[i] = c0 + c * (1.0 - exp(-distance/a));
            }
            break;
    }
    xVariogram[nx] = 1.0;

    return(xVariogram);
}

scalarList solveKriging
(
    const scalarList xList,
    const scalarList yList,
    const scalarList zList,
    const variogramType type,
    const scalar a,
    const scalar c,
    const scalar c0,
    const scalar x,
    const scalar y,
    const scalar z,
    SquareMatrix<scalar>& invRBFmatrix
)
{
    solveKriging(xList, yList, zList, type, a, c, c0, invRBFmatrix);

    const scalarList xVariogram(getVariogramVector(
        xList, yList, zList, type, a, c, c0, x, y, z
    ));

    return(invRBFmatrix * xVariogram);
}

scalar kriging
(
    const scalarList xList,
    const scalarList yList,
    const scalarList zList,
    const scalarList vList,
    const variogramType type,
    const scalar a,
    const scalar c,
    const scalar c0,
    const scalar x,
    const scalar y,
    const scalar z,
    SquareMatrix<scalar>& invRBFmatrix
)
{
    const scalarList xVariogram(getVariogramVector(
        xList, yList, zList, type, a, c, c0, x, y, z
    ));

    const scalarList w(invRBFmatrix * xVariogram);

    scalar res(0);
    forAll(vList, i)
    {
        res += w[i] * vList[i];
    }

    return(res);
}

scalar kriging
(
    const scalarList w,
    const scalarList vList
)
{
    scalar res(0);
    forAll(vList, i)
    {
        res += w[i] * vList[i];
    }

    return(res);
}

scalar kriging
(
    const scalarList xVariogram,
    const scalarList vList,
    SquareMatrix<scalar>& invRBFmatrix
)
{
    const scalarList w(invRBFmatrix * xVariogram);

    scalar res(0);
    forAll(vList, i)
    {
        res += w[i] * vList[i];
    }

    return(res);
}


}

}


// ************************************************************************* //
