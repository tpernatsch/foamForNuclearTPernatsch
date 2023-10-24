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

#include "latticeMap.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// namespace Foam
// {
//     defineTypeNameAndDebug(neutronics, 0);
//     defineRunTimeSelectionTable(neutronics, dictionary);
// }

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::latticeMap::latticeMap
(
    const dictionary& dict
)
:
    dict_(dict),
    type_(dict.get<word>("type")),
    originX_(dict.get<scalarList>("origin")[0]),
    originY_(dict.get<scalarList>("origin")[1]),
    pitch_(dict.get<scalar>("pitch")),
    rotationAngle_(dict.lookupOrDefault<scalar>("rotationAngle", 0)),
    nx_(dict.get<labelList>("nElements")[0]),
    ny_(dict.get<labelList>("nElements")[1]),
    lattice_(dict.get<wordList>("lattice")),
    xPos_(0),
    yPos_(0)
{
    if (nx_ * ny_ != lattice_.size())
    {
        FatalErrorInFunction
            << "nElements and lattice size are different"
            << exit(FatalError);
    }

    const scalar summitToSummit(pitch_ * 2.0/sqrt(3.0));
    const word nullElement("e");
    for (int i = 0; i < nx_; i++)
    {
        for (int j = 0; j < ny_; j++)
        {
            if (lattice_[ny_*i+j] != nullElement)
            {
                scalar x(0), y(0);
                if (type_ == "square")
                {
                    x = pitch_ * (i - (nx_-1.0)/2.0);
                    y = pitch_ * (j - (ny_-1.0)/2.0);
                }
                else if (type_ == "hexagonal")
                {
                    x = 0.75 * summitToSummit * (i - (nx_-1.0)/2.0);
                    y = pitch_ * ((j - (ny_-1.0)/2.0) + 0.5 * (i - (nx_-1.0)/2.0));
                }
                else
                {
                    FatalErrorInFunction
                        << "Lattice type " << type_ << " is not correct. Available lattice types: "
                        << "hexagonal or square."
                        << exit(FatalError);
                }
                xPos_.append(
                    originX_ + x*cos(rotationAngle_) - y*sin(rotationAngle_)
                );
                yPos_.append(
                    originY_ + x*sin(rotationAngle_) + y*cos(rotationAngle_)
                );
            }
        }
    }

    Info<< xPos_ << endl;
    Info<< yPos_ << endl;
}

Foam::latticeMap::latticeMap
(
    IOdictionary object,
    word latticeMapName
)
// :
//     dict_
//     type_
//     origin_
//     pitch_
//     nElements_
//     lattice_
    // startTime_(0.0),
    // functionPtr_(nullptr)
{
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


Foam::scalarList Foam::latticeMap::getXpositions() const
{
    return(xPos_);
}

Foam::scalarList Foam::latticeMap::getYpositions() const
{
    return(yPos_);
}


// ************************************************************************* //
