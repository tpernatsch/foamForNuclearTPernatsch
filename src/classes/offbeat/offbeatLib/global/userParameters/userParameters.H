/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Global
    userParameters

Description
    Dictionary reading and supplying user-modifiable parameters.

    The values are read from system/userParameters.

SourceFiles
    userParameters.C

\*---------------------------------------------------------------------------*/

#ifndef userParameters_H
#define userParameters_H

#include "dictionary.H"
#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Dictionary to record all registered user parameters
HashTable<dimensionedScalar, word>& userParametersRegistry();

//- Return a reference to the user parameters dictionary
HashTable<dimensionedScalar, word>& userParameters();

//- Register and read a user parameter or return a default value
dimensionedScalar userParameter
(
    const char* group,
    const char* name,
    const dimensionSet& dimensions,
    const scalar defaultValue
);

//- Return a user parameter
dimensionedScalar userParameter
(
    const char* group,
    const char* name
);

//- Return true if parameter is provided by the user
bool isInUserParameters
(
    const char* group,
    const char* name
);

//- Return the default value of a user parameter from the registry
dimensionedScalar registeredUserParameter
(
    const char* group,
    const char* name
);

//- Display a list of user-specified parameters
void listUserParameters();

//- Display a list of registered user parameters
void listRegisteredUserParameters();

//- Check the supplied user parameters in system/userParameters
void checkUserParameters();

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define defineParameter(classname, Type, name, dimensions, defaultValue)   \
const Foam::dimensionedScalar classname::Type                                   \
(                                                                               \
    Foam::userParameter                                                         \
    (                                                                           \
        classname::group_,                                                      \
        name,                                                                   \
        Foam::dimensionSet dimensions,                                          \
        defaultValue                                                            \
    )                                                                           \
);

# define declareParameterWithVirtual(classname, name, func)   \
static const dimensionedScalar name;                            \
virtual const dimensionedScalar& func() const                   \
{                                                               \
    return classname::name;                                     \
}   

# define declareNameWithVirtual(classname, name, func)   \
static const char* name;                                             \
virtual const char*& func() const                                    \
{                                                                   \
    return classname::name;                                         \
}                                                                              

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
