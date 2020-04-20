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

#include "phaseConstant.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(phaseConstant, 0);
    addToRunTimeSelectionTable(dragModel, phaseConstant, FSDragModels);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::phaseConstant::phaseConstant
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
    alpha_(FSPair.fluidRef()),
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

void Foam::dragModels::phaseConstant::correctKd(volTensorField& Kd) const
{   
    forAll(cellList_, i)
    {
        label celli(cellList_[i]);
        tensor value(alpha_[celli]*Kd_);
        
        tensor& Kdi(Kd[celli]);
        Kdi[0] = value[0];
        Kdi[4] = value[4];
        Kdi[8] = value[8];
    }
}


// ************************************************************************* //
