/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2312                                                  |
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
    Foam::nuclearDataOneEnergy

Description
    Subclass for nuclear data handling. The interpolation scheme used is based
    on the polyharmonic spline radial basis function RBF. The RBF used is
    linear.

\mainauthor of this file (and associated .C or included .H files):
    Thomas Guilbaud <thomas.guilbaud@epfl.ch>, EPFL/Transmutex SA (Switzerland)

SourceFiles
    nuclearDataOneEnergy.C

\*---------------------------------------------------------------------------*/

#ifndef nuclearDataOneEnergy_H
#define nuclearDataOneEnergy_H

#include "volFields.H"
#include "LUscalarMatrix.H"
#include "fvCFD.H"
#include "radialBasisFunctionInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class nuclearDataOneEnergy Declaration
\*---------------------------------------------------------------------------*/

class nuclearDataOneEnergy
{

private:

    //- Cellzone name
    word zoneName_;

    //- Number of perturbed state
    label nParameters_;

    //- Number of perturbed state found and valid
    // label nParametersReduced_;

    //-
    List<label> reducedParamIdx_;

    //- List of found perturbed state in nuclearData
    List<bool> isParameterList_;

    //- XS values
    scalarList data_;

    //- Parameter list corespond to one XS value
    List<scalarList> parameterList_;

    //- Perturbation data as a N data times N used parameters list
    List<scalarList> pointList_;

    //- Polyharmonique Spline Weights
    scalarList weights_;

    //- Inversed RBF matrix
    SquareMatrix<scalar> invRBFmatrix_;

    //-
    scalarList currValue_;

    //-
    scalarList currValueReduced_;

public:

    //- Constructors

        //- Constructor from dictionary
        nuclearDataOneEnergy(const word zoneName);

    //- Destructor
        virtual ~nuclearDataOneEnergy(){}


    //- Member functions

        //- Add data point for interpolation
        void addData(
            const scalar value,
            const scalar Tfuel,
            const scalar Tclad,
            const scalar Tcool,
            const scalar TstructMech,
            const scalar rhoCool,
            const scalar axExp,
            const scalar radExp
        );

        //- Build interpolation scheme
        void build();

        //-
        scalar getRef() const;

        //- Get the interpolated nuclear data value based on the input
        //  parameters
        scalar get(
            const scalar Tfuel,
            const scalar Tclad,
            const scalar Tcool,
            const scalar TstructMech,
            const scalar rhoCool,
            const scalar axExp,
            const scalar radExp,
            const bool isParametrize = true
        );

private:

    //- Private member functions

        //- Test if scalarList all same value, if true the interpolation on this
        //  variables is useless
        bool isAllSameValues(const scalarList list) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // nuclearDataOneEnergy

// ************************************************************************* //