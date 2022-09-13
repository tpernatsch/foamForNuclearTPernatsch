/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "timeProfile.H"

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
    const dictionary& dict
)
:
    dict_(dict),
    type_(dict_.get<word>("type")),
    startTime_(dict_.lookupOrDefault<scalar>("startTime", 0.0))
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

Foam::timeProfile::timeProfile
(
    IOdictionary object,
    word timeProfileName
)
:
    startTime_(0.0),
    functionPtr_(nullptr)
{
    if (object.found(timeProfileName))
    {
        dict_ = object.subDict(timeProfileName);

        type_ = dict_.get<word>("type");

        startTime_ = dict_.lookupOrDefault<scalar>("startTime", 0.0);

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
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::timeProfile::valid()
{
    return(functionPtr_.valid());
}

scalar Foam::timeProfile::value(scalar time)
{
    return
    (
        functionPtr_.valid()
            ? functionPtr_->value(time - startTime_)
            : 0.0
    );
}

// ************************************************************************* //
