/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2212                                                  |
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

    //- Reference values
    std::vector<scalar> xRef{0.1, 1.0, 2.0, 3.0};
    std::vector<scalar> yRef{1, 10.0, 20.0, 30.0};

    //- Loop over the tests
    int nFailed(0);
    int nCorrect(0);
    for (size_t i = 0; i < xRef.size(); i++)
    {
        scalarList x_i(1);
        x_i[0] = xRef[i];
        const scalar y_i(yRef[i]);
        const scalar y_e(polyharmonicSpline(weights, posList, x_i));

        //- Relative error test
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

    //- Resume all tests
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

    //- Reference values
    std::vector<scalar> xRef{ 1.0,  1.0,  2.0,  1.0,  2.0};
    std::vector<scalar> yRef{ 2.0,  4.0,  2.0,  8.0,  4.0};
    std::vector<scalar> zRef{10.0, 20.0, 40.0, 40.0, 50.0};

    //- Loop over the tests
    int nFailed(0);
    int nCorrect(0);
    for (size_t i = 0; i < xRef.size(); i++)
    {
        scalarList xy_i(2);
        xy_i[0] = xRef[i];
        xy_i[1] = yRef[i];
        const scalar z_i(zRef[i]);
        const scalar z_e(polyharmonicSpline(weights, posList, xy_i));

        //- Relative error test
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

    //- Resume all tests
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

    return 0;
}

// ************************************************************************* //
