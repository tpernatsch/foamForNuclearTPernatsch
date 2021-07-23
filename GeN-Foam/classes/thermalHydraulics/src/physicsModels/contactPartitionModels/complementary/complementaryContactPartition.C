/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "FSPair.H"
#include "complementaryContactPartition.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace contactPartitionModels
{
    defineTypeNameAndDebug(complementary, 0);
    addToRunTimeSelectionTable
    (
        contactPartitionModel,
        complementary,
        contactPartitionModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::contactPartitionModels::complementary::complementary
(
    const FSPair& pair,
    const dictionary& dict,
    const objectRegistry& objReg
)
:
    contactPartitionModel
    (
        pair,
        dict,
        objReg
    ),
    //complementaryPair_(nullptr),
    complementaryModel_(nullptr)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::contactPartitionModels::complementary::value
(
    const label& celli
) const
{
    if (complementaryModel_ == nullptr)
    {
        HashTable<const contactPartitionModel*> models
        (
            mesh_.lookupClass<contactPartitionModel>()
        );
        complementaryModel_ = 
        (
            (this->type() == models[models.toc()[0]]->type()) ?
            models[models.toc()[1]] : models[models.toc()[0]]
        );
    }

    return 1.0 - complementaryModel_->value(celli);
}

void Foam::contactPartitionModels::complementary::correctField
(
    volScalarField& f
)
{
    /*
    if (complementaryPair_ == nullptr)
    {
        HashTable<const FSPair*> pairs(mesh_.lookupClass<FSPair>());
        complementaryPair_ = 
        &(
            (pair_.name() == pairs[pairs.toc()[0]]->name()) ? 
            pairs[pairs.toc()[1]] : pairs[pairs.toc()[0]]
        );
    }
    f = 1.0-complementaryPair_->f();
    */
    if (complementaryModel_ == nullptr)
    {
        HashTable<const contactPartitionModel*> models
        (
            mesh_.lookupClass<contactPartitionModel>()
        );
        complementaryModel_ = 
        (
            (this->type() == models[models.toc()[0]]->type()) ?
            models[models.toc()[1]] : models[models.toc()[0]]
        );
    }

    forAll(mesh_.cells(), i)
    {
        f[i] = 1.0 - complementaryModel_->value(i);
    }
    f.correctBoundaryConditions();
}

// ************************************************************************* //
