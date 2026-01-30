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

#include "OTstream.H"
#include "radialBasisFunctionInterpolation.H"
#include "IOmanip.H"

using namespace Foam;
using namespace radialBasisFunctionInterpolation;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


scalar relativeError(scalar value, scalar th)
{
    return(2.0*mag((th-value)/(th+value)));
}


void printSummary(const int nFailed, const int nCorrect)
{
    if (nFailed == 0)
    {
        Info << "\033[32m" << "passed" << "\033[0m"
             << " with " << nCorrect << " points"
             << endl << endl;
    }
    else
    {
        Info << endl;
        Info << "  => \033[33mWARNING\033[0m: "
             << nFailed << "/" << nCorrect+nFailed << " tests failed"
             << endl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int test_PHSRBFlin1D(const scalar threshold)
{
    Info<< "Test 1D polyharmonic spline RBF with " << threshold*100 << "% threshold ... ";

    List<scalarList> posList(1);
    scalarList xList(2);
    xList[0] = 1.0;
    xList[1] = 2.0;
    posList[0] = xList;

    scalarList values(2);
    values[0] = 10.0;
    values[1] = 20.0;

    SquareMatrix<scalar> invRBFmatrix(0);

    scalarList weights = solvePolyharmonicSpline(posList, values, invRBFmatrix);

    //  Reference values
    std::vector<scalar> xRef{0.1, 1.0, 2.0, 3.0};
    std::vector<scalar> yRef{1, 10.0, 20.0, 30.0};

    //  Loop over the tests
    int nFailed(0);
    int nCorrect(0);
    for (size_t i = 0; i < xRef.size(); i++)
    {
        scalarList x_i(1);
        x_i[0] = xRef[i];
        const scalar y_i(yRef[i]);
        const scalar y_e(polyharmonicSpline(weights, posList, x_i));

        //  Relative error test
        const scalar relErr(relativeError(y_i, y_e));
        if (relErr > threshold)
        {
            Info<< endl;
            Info<< "    **** \033[33mWARNING\033[0m: Foam::polyharmonicSpline("
                << weights << ", "
                << posList << ", "
                << x_i
                << ") = " << y_e << " (exp: " << y_i
                << ") -> rel error = " << relErr;
            nFailed++;
        }
        else
        {
            nCorrect++;
        }
    }

    //  Resume all tests
    printSummary(nFailed, nCorrect);

    return nFailed;
}


int test_PHSRBFlin2D(const scalar threshold)
{
    Info<< "Test 2D polyharmonic spline RBF with " << threshold*100 << "% threshold ... ";

    List<scalarList> posList(2);
    scalarList xList(3);
    scalarList yList(3);
    xList[0] = 1.0; yList[0] = 2.0;
    xList[1] = 1.0; yList[1] = 4.0;
    xList[2] = 2.0; yList[2] = 2.0;
    posList[0] = xList;
    posList[1] = yList;

    scalarList values(3);
    values[0] = 10.0;
    values[1] = 20.0;
    values[2] = 40.0;

    SquareMatrix<scalar> invRBFmatrix(0);

    scalarList weights = solvePolyharmonicSpline(posList, values, invRBFmatrix);

    //  Reference values
    std::vector<scalar> xRef{ 1.0,  1.0,  2.0,  1.0,  2.0};
    std::vector<scalar> yRef{ 2.0,  4.0,  2.0,  8.0,  4.0};
    std::vector<scalar> zRef{10.0, 20.0, 40.0, 40.0, 50.0};

    //  Loop over the tests
    int nFailed(0);
    int nCorrect(0);
    for (size_t i = 0; i < xRef.size(); i++)
    {
        scalarList xy_i(2);
        xy_i[0] = xRef[i];
        xy_i[1] = yRef[i];
        const scalar z_i(zRef[i]);
        const scalar z_e(polyharmonicSpline(weights, posList, xy_i));

        //  Relative error test
        const scalar relErr(relativeError(z_i, z_e));
        if (relErr > threshold)
        {
            Info<< endl;
            Info<< "    **** \033[33mWARNING\033[0m: Foam::polyharmonicSpline("
                << weights << ", "
                << posList << ", "
                << xy_i
                << ") = " << z_e << " (exp: " << z_i
                << ") -> rel error = " << relErr;
            nFailed++;
        }
        else
        {
            nCorrect++;
        }
    }

    //  Resume all tests
    printSummary(nFailed, nCorrect);

    return nFailed;
}

int test_PHSRBFlin2DSmallValues(const scalar threshold)
{
    Info<< "Test 2D polyharmonic spline RBF for small values with " << threshold*100 << "% threshold ... ";

    List<scalarList> posList(2);
    const label nPoints(36);
    scalarList xList(nPoints);
    scalarList yList(nPoints);
    scalarList values(nPoints);
    xList[0] = 371.15; yList[0] = 989.16; values[0] = 9.019256e-12;

    xList[1] = 300.0; yList[1] = 600.0; values[1] = 8.745869e-12;
    xList[2] = 300.0; yList[2] = 700.0; values[2] = 8.844590e-12;
    xList[3] = 300.0; yList[3] = 800.0; values[3] = 8.916762e-12;
    xList[4] = 300.0; yList[4] = 900.0; values[4] = 8.977144e-12;
    xList[5] = 300.0; yList[5] = 1000.0; values[5] = 9.025965e-12;

    xList[6] = 400.0; yList[6] = 600.0; values[6] = 8.748343e-12;
    xList[7] = 400.0; yList[7] = 700.0; values[7] = 8.844299e-12;
    xList[8] = 400.0; yList[8] = 800.0; values[8] = 8.920262e-12;
    xList[9] = 400.0; yList[9] = 900.0; values[9] = 8.979478e-12;
    xList[10] = 400.0; yList[10] = 1000.0; values[10] = 9.023146e-12;

    xList[11] = 500.0; yList[11] = 600.0; values[11] = 8.748084e-12;
    xList[12] = 500.0; yList[12] = 700.0; values[12] = 8.842483e-12;
    xList[13] = 500.0; yList[13] = 800.0; values[13] = 8.921235e-12;
    xList[14] = 500.0; yList[14] = 900.0; values[14] = 8.977545e-12;
    xList[15] = 500.0; yList[15] = 1000.0; values[15] = 9.028591e-12;

    xList[16] = 600.0; yList[16] = 600.0; values[16] = 8.746593e-12;
    xList[17] = 600.0; yList[17] = 700.0; values[17] = 8.843715e-12;
    xList[18] = 600.0; yList[18] = 800.0; values[18] = 8.920047e-12;
    xList[19] = 600.0; yList[19] = 900.0; values[19] = 8.975632e-12;
    xList[20] = 600.0; yList[20] = 1000.0; values[20] = 9.024377e-12;

    xList[21] = 700.0; yList[21] = 600.0; values[21] = 8.745491e-12;
    xList[22] = 700.0; yList[22] = 700.0; values[22] = 8.839599e-12;
    xList[23] = 700.0; yList[23] = 800.0; values[23] = 8.916838e-12;
    xList[24] = 700.0; yList[24] = 900.0; values[24] = 8.975146e-12;
    xList[25] = 700.0; yList[25] = 1000.0; values[25] = 9.028591e-12;

    xList[26] = 800.0; yList[26] = 600.0; values[26] = 8.743773e-12;
    xList[27] = 800.0; yList[27] = 700.0; values[27] = 8.837589e-12;
    xList[28] = 800.0; yList[28] = 800.0; values[28] = 8.921473e-12;
    xList[29] = 800.0; yList[29] = 900.0; values[29] = 8.977415e-12;
    xList[30] = 800.0; yList[30] = 1000.0; values[30] = 9.026030e-12;

    xList[31] = 900.0; yList[31] = 600.0; values[31] = 8.741957e-12;
    xList[32] = 900.0; yList[32] = 700.0; values[32] = 8.843273e-12;
    xList[33] = 900.0; yList[33] = 800.0; values[33] = 8.916871e-12;
    xList[34] = 900.0; yList[34] = 900.0; values[34] = 8.975114e-12;
    xList[35] = 900.0; yList[35] = 1000.0; values[35] = 9.026678e-12;
    posList[0] = xList;
    posList[1] = yList;

    SquareMatrix<scalar> invRBFmatrix(0);

    scalarList weights = solvePolyharmonicSpline(posList, values, invRBFmatrix);

    //  Reference values
    std::vector<scalar> xRef{301.0, 350.0};
    std::vector<scalar> yRef{600.0, 950.0};
    std::vector<scalar> zRef{8.745962599638433e-12, 9.001298609220322e-12};

    //  Loop over the tests
    int nFailed(0);
    int nCorrect(0);
    for (size_t i = 0; i < xRef.size(); i++)
    {
        scalarList xy_i(2);
        xy_i[0] = xRef[i];
        xy_i[1] = yRef[i];
        const scalar z_i(zRef[i]);
        const scalar z_e(polyharmonicSpline(weights, posList, xy_i));

        //  Relative error test
        const scalar relErr(relativeError(z_i, z_e));
        if (relErr > threshold)
        {
            Info<< endl;
            Info<< "    **** \033[33mWARNING\033[0m: Foam::polyharmonicSpline("
                << weights << ", "
                << posList << ", "
                << xy_i
                << ") = " << z_e << " (exp: " << z_i
                << ") -> rel error = " << relErr;
            nFailed++;
        }
        else
        {
            nCorrect++;
        }
    }

    //  Resume all tests
    printSummary(nFailed, nCorrect);

    return nFailed;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    if (test_PHSRBFlin1D(0.0001) >= 1)
    {
        return 1;
    }
    if (test_PHSRBFlin2D(0.0001) >= 1)
    {
        return 1;
    }
    if (test_PHSRBFlin2DSmallValues(0.0002) >= 1)
    {
        return 1;
    }

    return 0;
}

// ************************************************************************* //
