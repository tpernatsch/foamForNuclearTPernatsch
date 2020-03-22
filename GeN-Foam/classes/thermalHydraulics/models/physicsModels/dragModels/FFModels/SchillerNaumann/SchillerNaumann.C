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

#include "SchillerNaumann.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(SchillerNaumann, 0);
    addToRunTimeSelectionTable(dragModel, SchillerNaumann, FFDragModels);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::SchillerNaumann::SchillerNaumann
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
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dragModels::SchillerNaumann::correctKd(volTensorField& Kd) const
{
    const volScalarField& Re(FFPair_->Re());
    const volScalarField& alpha1(FFPair_->fluid1());
    const volScalarField& alpha2(FFPair_->fluid2());
    const volScalarField& magUr(FFPair_->magUr());
    const volScalarField& DhDispersed(FFPair_->DhDispersed());
    const volScalarField& rhoContinuous(FFPair_->rhoContinuous());

    //- I don't care about mesh boundaries, cell-by-cell is faster
    forAll(mesh_.cells(), i)
    {
        tensor& Kdi(Kd[i]);
        scalar value
        (
            0.75*alpha1[i]*alpha2[i]*rhoContinuous[i]*magUr[i]/DhDispersed[i]
        );
        const scalar& Rei(Re[i]);
        if (Rei < 1000) value *= 24*(1.0+0.15*pow(Rei, 0.687))/Rei;
        else value *= 0.44;
        Kdi[0] = value;
        Kdi[4] = value;
        Kdi[8] = value;
    }
}


// ************************************************************************* //
