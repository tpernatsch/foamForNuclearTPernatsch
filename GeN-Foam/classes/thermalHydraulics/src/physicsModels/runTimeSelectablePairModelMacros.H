/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2015 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#ifndef runTimeSelectablePairModel_H
#define runTimeSelectablePairModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define STRING(s) #s

//- Declare class className for a pair model with pair of type pairType
#define declareRunTimeSelectablePairModel(pairType, className)              \
class className : IOdictionary                                              \
{                                                                           \
protected:                                                                  \
    const pairType& pair_;                                                  \
public:                                                                     \
    TypeName(STRING(className));                                            \
    declareRunTimeSelectionTable                                            \
    (                                                                       \
        autoPtr,                                                            \
        className,                                                          \
        className##s,                                                       \
        (                                                                   \
            const pairType& pair,                                           \
            const dictionary& dict,                                         \
            const objectRegistry& objReg                                    \
        ),                                                                  \
        (pair, dict, objReg)                                                \
    );                                                                      \
    className                                                               \
    (                                                                       \
        const pairType& pair,                                               \
        const dictionary& dict,                                             \
        const objectRegistry& objReg                                        \
    );                                                                      \
    virtual ~className(){}                                                  \
    static autoPtr<className> New                                           \
    (                                                                       \
        const pairType& pair,                                               \
        const dictionary& dict,                                             \
        const objectRegistry& objReg                                        \
    );                                                                      \
    virtual scalar value(const label& celli) const = 0;                     \
};

//- Define class className for a pair model with pair of type pairType
#define defineRunTimeSelectableBasePairModel(pairType, className)           \
namespace Foam                                                              \
{                                                                           \
    defineTypeNameAndDebug(className, 0);                                   \
    defineRunTimeSelectionTable                                             \
    (                                                                       \
        className,                                                          \
        className##s                                                        \
    );                                                                      \
}                                                                           \
Foam::autoPtr<className> Foam::className::New                               \
(                                                                           \
    const pairType& pair,                                                   \
    const dictionary& dict,                                                 \
    const objectRegistry& objReg                                            \
)                                                                           \
{                                                                           \
    word type(dict.get<word>("type"));                                      \
    Info<< "Selecting " << STRING(className) << " for pair "                \
        << pair.name() << ": " << type << endl;                             \
    className##sConstructorTable::iterator cstrIter =                       \
        className##sConstructorTablePtr_->find(type);                       \
    if (cstrIter == className##sConstructorTablePtr_->end())                \
    {                                                                       \
        FatalErrorInFunction                                                \
            << "Unknown " << STRING(className) << " type "                  \
            << type << endl << endl                                         \
            << "Valid " << STRING(className) << " types "                   \
            << "are: " << endl                                              \
            << className##sConstructorTablePtr_->sortedToc()                \
            << exit(FatalError);                                            \
    }                                                                       \
    return cstrIter()(pair, dict, objReg);                                  \
}                                                                           \
Foam::className::className                                                  \
(                                                                           \
    const pairType& pair,                                                   \
    const dictionary& dict,                                                 \
    const objectRegistry& objReg                                            \
)                                                                           \
:                                                                           \
    IOdictionary                                                            \
    (                                                                       \
        IOobject                                                            \
        (                                                                   \
            STRING(className),                                              \
            pair.mesh().time().timeName(),                                  \
            objReg,                                                         \
            IOobject::NO_READ,                                              \
            IOobject::NO_WRITE                                              \
        ),                                                                  \
        dict                                                                \
    ),                                                                      \
    pair_(pair)                                                             \
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
