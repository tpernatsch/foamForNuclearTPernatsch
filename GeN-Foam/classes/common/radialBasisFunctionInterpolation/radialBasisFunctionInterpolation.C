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

        // solve(w, A, vList);
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


scalarList solvePolyharmonicSpline
(
    const scalarList xList,
    const scalarList yList,
    const scalarList zList,
    const scalarList vList,
    SquareMatrix<scalar>& invRBFmatrix
)
{
    const label nx(xList.size());

    scalarList w(nx+4, 0.0);

    if (invRBFmatrix.m() != nx+4)
    {
        invRBFmatrix.resize(nx+4);

        SquareMatrix<scalar> A(nx+4, 0.0);

        forAll(xList, i)
        {
            forAll(xList, j)
            {
                if (i != j)
                {
                    const scalar r(sqrt(
                        sqr(xList[i]-xList[j]) 
                        + sqr(yList[i]-yList[j])
                        + sqr(zList[i]-zList[j])
                    ));
                    A[i][j] = sqr(r) * log(r);
                }
                else
                {
                    A[i][j] = 0;
                }
            }
            // Polynomial correction
            A[i][nx] = 1;
            A[i][nx+1] = xList[i];
            A[i][nx+2] = yList[i];
            A[i][nx+3] = zList[i];
            A[nx][i] = 1;
            A[nx+1][i] = xList[i];
            A[nx+2][i] = yList[i];
            A[nx+3][i] = zList[i];
        }

        // Inverse the matrix once and store it for later iterations
        LUscalarMatrix Atemp(A);

        Atemp.inv(invRBFmatrix);

        // solve(w, A, vListTemp);
    }

    scalarList vListTemp(vList);
    for (int i = 0; i < 4; i++)
    {
        vListTemp.append(0);
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
    const scalar z
)
{
    const label nx(xList.size());
    scalar res(0);
    forAll(xList, i)
    {
        const scalar r(sqrt(
            sqr(x-xList[i]) 
            + sqr(y-yList[i])
            + sqr(z-zList[i])
        ));
        res += w[i] * sqr(r) * log(r);
    }
    res += w[nx] + w[nx+1]*x + w[nx+2]*y + w[nx+3]*z;
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
            const scalar rd(sqrt(0*sqr(dx) + sqr(dy)));
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
                const scalar dphix((2.0 * dx * log(rd) + dx));
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
    forAll(xList, i)
    {
        const scalar r(sqrt(
            sqr(x-xList[i]) 
            + sqr(y-yList[i])
            + sqr(z-zList[i])
        ));
        if (r > 0)
        {
            res += (
                w[i] + w[i+nx] * (x-xList[i]) + w[i+2*nx] * (y-yList[i])
            ) * sqr(r) * log(r);
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
    SquareMatrix<scalar>& invRBFmatrix
)
{
    const label nx(xList.size());

    scalarList w(nx+1, 0.0);

    if (invRBFmatrix.m() != nx+1)
    {
        invRBFmatrix.resize(nx+1);

        SquareMatrix<scalar> A(nx+1, 0.0);

        // const labelList& regionCells(structure_.cellLists()[region]);
        const scalarList& V(mesh.V());

        forAll(xList, i)
        {
            forAll(xList, j)
            {
                if (i != j)
                {
                    const scalar r(sqrt(
                        sqr(xList[i]-xList[j]) 
                        + sqr(yList[i]-yList[j])
                        + sqr(zList[i]-zList[j])
                    ));
                    A[i][j] = sqr(r) * log(r);
                }
                else
                {
                    A[i][j] = 0;
                }
            }
            // Polynomial correction
            A[i][nx] = 1;
            // A[i][nx+1] = xList[i];
            // A[i][nx+2] = yList[i];
            // A[i][nx+3] = zList[i];
            // A[i][nx+1] = zList[i];
            // A[nx][i] = 1;
            // A[nx+1][i] = xList[i];
            // A[nx+2][i] = yList[i];
            // A[nx+3][i] = zList[i];
            // A[nx+1][i] = zList[i];

            scalar totalVolume(0), totalPhi(0);
            forAll(regionCells, celli)
            {
                const scalar xCell(mesh.C().internalField()[celli].x());
                const scalar yCell(mesh.C().internalField()[celli].y());
                const scalar zCell(mesh.C().internalField()[celli].z());
                const scalar xPos(xCell-xList[i]);
                const scalar yPos(yCell-yList[i]);
                const scalar zPos(zCell-zList[i]);
                // label cellNumber = mesh.findCell(point(xPos, yPos, zPos));
                // totalPower += //alpha_[celli]
                //     /***/ fractionOfPowerFromNeutronics_[regioni]
                //     * structure_.powerDensityNeutronics()[celli]
                //     * V[celli];
                totalVolume += V[celli];
                const scalar rsqr(sqr(xPos)+sqr(yPos)+sqr(zPos));
                // totalPhi += (rsqr * log(sqrt(rsqr))) * V[cellNumber];
                totalPhi += (rsqr * log(sqrt(rsqr))) * V[celli];
            }

            // Info<< totalPhi << " " << totalVolume << endl;

            A[nx][i] = totalPhi;
            A[nx][nx] = totalVolume;
        }

        // solve(w, A, vListTemp);

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
    const scalar z
)
{
    const label nx(xList.size());
    scalar res(0);
    forAll(xList, i)
    {
        const scalar r(sqrt(
            sqr(x-xList[i]) 
            + sqr(y-yList[i])
            + sqr(z-zList[i])
        ));
        res += w[i] * sqr(r) * log(r);
    }
    res += w[nx];
    return(res);
}



void solveKriging
(
    const scalarList xList,
    const scalarList yList,
    const scalarList zList,
    SquareMatrix<scalar>& invRBFmatrix
)
{
    const label nx(xList.size());

    if (invRBFmatrix.m() != nx+1)
    {
        invRBFmatrix.resize(nx+1);

        SquareMatrix<scalar> A(nx+1, 0.0);

        const scalar a(0.1), b(7.5), c(2.5);

        forAll(xList, i)
        {
            forAll(xList, j)
            {
                const scalar distance(sqrt
                (
                    sqr(xList[i]-xList[j]) + 
                    sqr(yList[i]-yList[j]) + 
                    sqr(zList[i]-zList[j])
                ));
                A[i][j] = c + b*(1.5*distance/a - 0.5*pow(distance/a, 3.0));
            }
            A[i][nx] = 1.0;
            A[nx][i] = 1.0;
        }

        LUscalarMatrix Atemp(A);

        Atemp.inv(invRBFmatrix);
    }
}


scalar kriging
(
    const scalarList xList,
    const scalarList yList,
    const scalarList zList,
    const scalarList vList,
    const scalar x,
    const scalar y,
    const scalar z,
    SquareMatrix<scalar>& invRBFmatrix
)
{
    const label nx(xList.size());

    const scalar a(0.1), b(7.5), c(2.5);
    scalarList xVariogram(nx+1);
    forAll(xList, i)
    {
        const scalar distance(sqrt
        (
            sqr(xList[i]-x) + 
            sqr(yList[i]-y) + 
            sqr(zList[i]-z)
        ));
        xVariogram[i] = c + b*(1.5*distance/a - 0.5*pow(distance/a, 3.0));
    }
    xVariogram[nx] = 1.0;

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
