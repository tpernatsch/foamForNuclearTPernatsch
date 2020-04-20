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

#include "constantDrag.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(constant, 0);
    addToRunTimeSelectionTable(dragModel, constant, FFDragModels);
    addToRunTimeSelectionTable(dragModel, constant, FSDragModels);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::constant::constant
(
    const objectRegistry& objReg,
    const dictionary& dict,
    const FFPair& FFPair
)
:
    dragModel
    (
        objReg,
        dict,
        FFPair
    ),
    Kd_
    (
        tensor
        (
            dict.get<scalar>("value"), 0, 0,
            0, dict.get<scalar>("value"), 0,
            0, 0, dict.get<scalar>("value")
        )
    )
{
    cellField_ = scalarField(mesh_.cells().size(), 0);
    forAll(mesh_.cells(), i)
    {
        cellList_.append(i);
        cellField_[i] = 1.0;
    }
}

Foam::dragModels::constant::constant
(
    const objectRegistry& objReg,
    const dictionary& dict,
    const FSPair& FSPair
)
:
    dragModel
    (
        objReg,
        dict,
        FSPair
    ),
    Kd_(tensor::zero)
{
    const entry* eKd(this->findEntry("value"));
    if (eKd->stream().size() == 1)
    {
        scalar Kd(this->get<scalar>("value"));
        Kd_[0] = Kd;
        Kd_[4] = Kd;
        Kd_[8] = Kd;
    }
    else
    {
        Kd_ = this->get<tensor>("value");
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dragModels::constant::correctKd(volTensorField& Kd) const
{   
    forAll(cellList_, i)
    {
        Kd[cellList_[i]] = Kd_;
    }
}


// ************************************************************************* //
