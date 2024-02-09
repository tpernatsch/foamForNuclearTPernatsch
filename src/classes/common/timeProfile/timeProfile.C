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

\*---------------------------------------------------------------------------*/

#if defined __has_include
#  if __has_include(<commDataLayer.H>) 
#    include <commDataLayer.H>
#    define isCommDataLayerIncluded
#  endif
#endif

#include "timeProfile.H"

#ifdef isCommDataLayerIncluded
#include "commDataLayer.H"
#endif

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// namespace Foam
// {
//     defineTypeNameAndDebug(neutronics, 0);
//     defineRunTimeSelectionTable(neutronics, dictionary);
// }

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeProfile::timeProfile
(
    const dictionary& dict,
    const Time& runTime
)
:
    dict_(dict),
    runTime_(runTime),
    type_(dict_.get<word>("type")),
    startTime_(dict_.lookupOrDefault<scalar>("startTime", 0.0)),
    functionPtr_(nullptr),
    isfmiPortSet_(false)
{
    if (type_ != "fmi")
    {
        functionPtr_.reset
        (
            Function1<scalar>::New
            (
                type_,
                dict_,
                type_
            )
        );
    }
    else
    {
        initializeFMI();
    }
}

Foam::timeProfile::timeProfile
(
    IOdictionary object,
    word timeProfileName,
    const Time& runTime
)
:
    runTime_(runTime),
    startTime_(0.0),
    functionPtr_(nullptr),
    isfmiPortSet_(false)
{
    if (object.found(timeProfileName))
    {
        dict_ = object.subDict(timeProfileName);
        
        type_ = dict_.get<word>("type");

        startTime_ = dict_.lookupOrDefault<scalar>("startTime", 0.0);

        if (type_ != "fmi")
        {
            functionPtr_.reset
            (
                Function1<scalar>::New
                (
                    type_,
                    dict_,
                    type_
                )
            );
        }
        else
        {
            initializeFMI();
        }
    }
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::timeProfile::initializeFMI()
{
#ifdef isCommDataLayerIncluded

    word keyFromFMU("nameFromFMU");
    if (dict_.found(keyFromFMU))
    {
        nameFromFMU_ = dict_.get<word>(keyFromFMU);

        Info<< dict_.dictName()
            << " FMI input name: " << nameFromFMU_ 
            << endl;

        // Communicating with the FMU
        commDataLayer& data = commDataLayer::New(runTime_); 
        // Store in data layer and set its initial value to 0 or user defined
        data.storeObj(
            dict_.lookupOrDefault<scalar>("initialValue", 0.0),
            nameFromFMU_,
            commDataLayer::causality::in
        );

        isfmiPortSet_ = true;
    }

#else // !isCommDataLayerIncluded

    FatalErrorInFunction
        << "FMI mode selected but FMI lib not built or linked to GeN-Foam."
        << exit(FatalError);

#endif // isCommDataLayerIncluded
}

bool Foam::timeProfile::valid() const
{
    if (isfmiPortSet_)
    {
        return(true);
    }
    if (!functionPtr_) // If functionPtr_ is not defined return false
    {
        return(false);
    }
    return(functionPtr_.valid());
}

scalar Foam::timeProfile::value(scalar time) const
{
#ifdef isCommDataLayerIncluded
    if (isfmiPortSet_)
    {
        commDataLayer& data = commDataLayer::New(runTime_);
        scalar value = data.getObj<scalar>(
            nameFromFMU_,
            commDataLayer::causality::in
        );
        return(value);
    }
#endif // isCommDataLayerIncluded

    return
    (
        functionPtr_.valid()
            ? functionPtr_->value(time - startTime_)
            : 0.0
    );
}

// ************************************************************************* //
