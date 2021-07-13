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

#include "twoPhaseDragMultiplierModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(twoPhaseDragMultiplierModel, 0);
    defineRunTimeSelectionTable
    (
        twoPhaseDragMultiplierModel, 
        twoPhaseDragMultiplierModels
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseDragMultiplierModel::twoPhaseDragMultiplierModel
(
    const fvMesh& mesh,
    const dictionary& dict,
    const objectRegistry& objReg
)
:
    IOdictionary
    (
        IOobject
        (
            typeName,
            mesh.time().timeName(),
            objReg,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dict
    ),
    mesh_(mesh),
    mFluidPtr_(nullptr),
    oFluidPtr_(nullptr),
    mKd_(nullptr),
    maxPhi2_
    (
        dict.lookupOrDefault<scalar>("maxPhi2", 1000)
    )
{
    //- Set ptr to multiplier fluid
    mFluidPtr_ = 
        &(
            mesh_.lookupObject<fluid>
            (
                word("alpha."+this->get<word>("multiplierFluid"))
            )
        );

    /*-----------------------------------------------------------------------*\
    |                                                                         |
    |                          H I C   M A N E B I T                          |  
    \*                                                                       */
    /*-----------------------------------------------------------------------*\
    |                                                                         |
    |   September 2020 was such a shit month it deserves a special mention.   |                                         
    |   In particular, Monday, September 7th 2020 saw a pretty damn bad bug   |
    |   down here which broke like 3 months of calculations without me        |
    |   realizing. Let us not forget Friday, September 25th 2020, which was   |
    |   just as shit, just not due to the computer code but due to, well,     |
    |   a certain Svetlana in my (now past) life. Sic transit!                |                                                          
    |                                                                         |
    |    - virmodoetiae                                                       |
    |                                                                         |
    \*-----------------------------------------------------------------------*/

    //- PORCO DIOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    
    //- Set ptr to other fluid
    HashTable<const fluid*> fluids(mesh_.lookupClass<fluid>());
    wordList fluidNames(fluids.toc());
    oFluidPtr_ = (fluidNames[0] == "alpha."+mFluidPtr_->name()) ?
        fluids[fluidNames[1]] : fluids[fluidNames[0]];

    //- Set ptr to mKd
    const FSPair& pair
    (
        mesh_.lookupObject<FSPair>(mFluidPtr_->name()+".structure")
    );
    mKd_ = &(pair.Kd());
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::twoPhaseDragMultiplierModel::onePhase(const label& celli) const
{
    return (mFluidPtr_->normalized()[celli] > 0.9999);
}

//- Default value
Foam::scalar Foam::twoPhaseDragMultiplierModel::phi2
(
    const label& celli
) const
{
    return 1.0;
}

Foam::tensor Foam::twoPhaseDragMultiplierModel::value
(
    const label& celli
) const
{
    return
        tensor
        (
            min(maxPhi2_, phi2(celli))*((*mKd_)[celli])*
            sqr(max(mFluidPtr_->normalized()[celli], 1e-9))*
            mFluidPtr_->magU()[celli]
        );
}

void Foam::twoPhaseDragMultiplierModel::correctField(volTensorField& KdTotU)
{
    forAll(mesh_.cells(), i)
    {
        KdTotU[i] = this->value(i);
    }
    KdTotU.correctBoundaryConditions();
}

// ************************************************************************* //
