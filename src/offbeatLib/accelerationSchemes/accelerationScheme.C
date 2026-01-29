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

#include "accelerationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
accelerationScheme<Type>::accelerationScheme
(
    label size,
    const dictionary& dict
)
:
    psi_(size),
    nFields_(0)
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class Type>
autoPtr<accelerationScheme<Type> > accelerationScheme<Type>::New
(
    label size,
    const dictionary& dict
)
{
    word schemeType(dict.lookup("type"));
        
    auto cstrIter
    = dictionaryConstructorTablePtr_->find(schemeType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn(__FUNCTION__)
            << "Unknown solution acceleration type " << schemeType
            << endl << endl
            << "Valid types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    Info<< "Selecting solution acceleration method "
        << schemeType << endl;

    return cstrIter()(size, dict);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
accelerationScheme<Type>::~accelerationScheme()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


template<class Type>
void accelerationScheme<Type>::addField(volTypeField& psi)
{
    if (nFields_ == psi_.size())
    {
        FatalErrorIn("accelerationScheme::addField")
            << "Attempt to add more than " << psi_.size() << "fields to "
            << "accelerationScheme." << abort(FatalError);
    }

    psi_.set(nFields_, &psi);
    nFields_++;
}


template<class Type>
void accelerationScheme<Type>::reset()
{}


template<class Type>
void accelerationScheme<Type>::accelerate()
{
    if (nFields_ != psi_.size())
    {
        FatalErrorIn("accelerationScheme::accelerate")
            << "Not all fields for acceleration are defined."
            << abort(FatalError);
    }
}


template<class Type>
int accelerationScheme<Type>::minIter() const
{
    return 1;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
