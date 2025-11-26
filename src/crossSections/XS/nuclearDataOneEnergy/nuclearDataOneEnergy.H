/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2506                                                  |
|    Copyright 2011-2016 OpenFOAM Foundation, 2017-2025 OpenCFD Ltd.         |
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
    Object for a single nuclear quantity handling at a specific energy. The
    interpolation scheme used is based on the polyharmonic spline radial basis
    function (RBF). This class is able to interpolate any
    number of perturbed states for any number of parameters. This class is mainly
    used by the ``XS`` object.

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

    //- Polyharmonic spline mode
    label polyharmonicSplineMode_;

    //- Number of perturbed state
    label nParameters_;

    //- List of indices of used perturb parameters
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

    //- Buffer list containing the input parameters
    // scalarList currValue_;

    //- Reduced buffer list containing the input parameters, only contains the
    //  perturbed parameters
    scalarList currValueReduced_;

public:

    //- Constructors

        //- Constructor from dictionary
        nuclearDataOneEnergy
        (
            const word zoneName,
            const label mode = 1,
            const label nParameters = 1
        );

    //- Destructor
        virtual ~nuclearDataOneEnergy(){}


    //- Member functions

        //- Add data point for interpolation
        void addData
        (
            const scalar value,
            const List<scalar>& parameters
        );

        //- Build interpolation scheme
        void build();

        //- Return the reference XS value
        scalar getRef() const;

        //- Get the interpolated nuclear data value based on the input
        //  parameters
        scalar get
        (
            const List<scalar>& parameters,
            const bool isParametrize = true
        );

private:

    //- Private member functions

        //- Test if scalarList all same value, if true the interpolation on this
        //- variables is useless
        bool isAllSameValues(const scalarList list) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // nuclearDataOneEnergy

// ************************************************************************* //