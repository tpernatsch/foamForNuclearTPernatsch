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

\*---------------------------------------------------------------------------*/


#include "nuclearDataOneEnergy.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nuclearDataOneEnergy::nuclearDataOneEnergy
(
    const word zoneName,
    const label mode,
    const label nParameters
) :
    zoneName_(zoneName),
    polyharmonicSplineMode_(mode),
    nParameters_(nParameters/*7*/),
    reducedParamIdx_(0),
    isParameterList_(nParameters_),
    data_(0),
    parameterList_(nParameters_),
    pointList_(0),
    weights_(0),
    // currValue_(nParameters_),
    currValueReduced_(0)
{

}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::nuclearDataOneEnergy::addData
(
    const scalar value,
    const List<scalar>& parameters
    // const scalar Tfuel,
    // const scalar Tclad,
    // const scalar Tcool,
    // const scalar TstructMech,
    // const scalar rhoCool,
    // const scalar axExp,
    // const scalar radExp
)
{
    data_.append(value);
    forAll(parameters, paramI)
    {
        parameterList_[paramI].append(parameters[paramI]);
    }
    // parameterList_[0].append(Tfuel);
    // parameterList_[1].append(Tclad);
    // parameterList_[2].append(Tcool);
    // parameterList_[3].append(TstructMech);
    // parameterList_[4].append(rhoCool);
    // parameterList_[5].append(axExp);
    // parameterList_[6].append(radExp);
}


void Foam::nuclearDataOneEnergy::build()
{
    // Look for valid perturbed parameters
    reducedParamIdx_.clear();
    forAll(parameterList_, paramI)
    {
        isParameterList_[paramI] = !isAllSameValues(parameterList_[paramI]);
        if (isParameterList_[paramI])
        {
            reducedParamIdx_.append(paramI);
        }
    }

    currValueReduced_.resize(reducedParamIdx_.size());

    // Build interpolation scheme
    if (data_.size() >= 2 && reducedParamIdx_.size() >= 1)
    {
        // Create reduced data point space, otherwise RBF matrix is non
        // inversible
        pointList_.clear();
        forAll(isParameterList_, paramI)
        {
            if (isParameterList_[paramI])
            {
                pointList_.append(parameterList_[paramI]);
            }
        }

        // Compute radial basis function weights
        weights_ = Foam::radialBasisFunctionInterpolation::solvePolyharmonicSpline
        (
            pointList_, data_, invRBFmatrix_, polyharmonicSplineMode_
        );
    }
}

scalar Foam::nuclearDataOneEnergy::getRef() const
{
    return(data_.first());
}

scalar Foam::nuclearDataOneEnergy::get
(
    // const scalar Tfuel,
    // const scalar Tclad,
    // const scalar Tcool,
    // const scalar TstructMech,
    // const scalar rhoCool,
    // const scalar axExp,
    // const scalar radExp,
    const List<scalar>& parameters,
    const bool isParametrize
)
{
    if (reducedParamIdx_.size() == 0 || !isParametrize)
    {
        return(data_.first());
    }

    // Build full parameter list
    // currValue_[0] = Tfuel;
    // currValue_[1] = Tclad;
    // currValue_[2] = Tcool;
    // currValue_[3] = TstructMech;
    // currValue_[4] = rhoCool;
    // currValue_[5] = axExp;
    // currValue_[6] = radExp;

    // Build reduced parameter list
    forAll(reducedParamIdx_, paramI)
    {
        // currValueReduced_[paramI] = currValue_[reducedParamIdx_[paramI]];
        currValueReduced_[paramI] = parameters[reducedParamIdx_[paramI]];
    }

    return
    (
        Foam::radialBasisFunctionInterpolation::polyharmonicSpline
        (
            weights_, pointList_, currValueReduced_, polyharmonicSplineMode_
            // weights_, pointList_, parameters, polyharmonicSplineMode_
        )
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::nuclearDataOneEnergy::isAllSameValues(const scalarList list) const
{
    const scalar first(list.first());
    forAll(list, idxI)
    {
        if (list[idxI] != first)
        {
            return(false);
        }
    }
    return(true);
}


// ************************************************************************* //
