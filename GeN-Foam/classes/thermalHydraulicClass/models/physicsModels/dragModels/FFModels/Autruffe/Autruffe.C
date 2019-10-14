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

#include "Autruffe.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(Autruffe, 0);
    addToRunTimeSelectionTable(dragModel, Autruffe, FFDragModels);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::Autruffe::Autruffe
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
    FFPair_
    (
        FFPair
    ),
    vapour_
    (
        mesh_.lookupObject<fluid>
        (
                "alpha."
            +   word(this->lookup("vapourName"))
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volTensorField> Foam::dragModels::Autruffe::Kd() const
{
    volScalarField alpha(vapour_/FFPair_.alphaSum());

    volScalarField Kds
    (
        4.31/(2*FFPair_.DhContinuous())*
        vapour_.rho()*
        pow
        (
            (1.0-alpha)*
            (1.0+75.0*(1.0-alpha)),
            0.95
        )
    );

    return 
        tmp<volTensorField>
        (
            new volTensorField(tensor::I*Kds)
        );
}


// ************************************************************************* //
