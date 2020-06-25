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

#include "NoKazimiInterfacialArea.H"
#include "fvCFD.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fluidInterfacialAreaModels
{
    defineTypeNameAndDebug(NoKazimiInterfacialArea, 0);
    addToRunTimeSelectionTable
    (
        fluidInterfacialAreaModel, 
        NoKazimiInterfacialArea, 
        fluidInterfacialAreaModels);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidInterfacialAreaModels::NoKazimiInterfacialArea::NoKazimiInterfacialArea
(
    const objectRegistry& objReg,
    const dictionary& dict,
    const fluid& dispersed,
    const fluid& continuous
)
:
    fluidInterfacialAreaModel
    (
        objReg,
        dict,
        dispersed,
        continuous
    ),
    alpha_
    (
        (dispersed.isGas()) ? dispersed : continuous
    ),
    D_("pinDiameter", dimLength, dict),
    P_("pinPitch", dimLength, dict),
    PD_(P_/D_),
    A_
    (
        4.0*constant::mathematical::pi/
        (
            D_*
            (2*Foam::sqrt(3.0)*sqr(PD_)-constant::mathematical::pi)
        )
    )
{
    if 
    (
        dispersed.isLiquid() and continuous.isLiquid()
    )
    {
        FatalErrorInFunction
            << "The stateOfMatter of one of the two phases must be gas!"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidInterfacialAreaModels::NoKazimiInterfacialArea::~NoKazimiInterfacialArea()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> 
Foam::fluidInterfacialAreaModels::NoKazimiInterfacialArea::iA() const
{
    volScalarField alphaSum(dispersed_+continuous_);
    
    volScalarField alpha(alpha_/alphaSum);

    tmp<volScalarField> tiA
    (
        new volScalarField
        (
            A_*alpha*
            min
            (
                (1.0-alpha)/0.043,
                1.0
            )
        )
    );
    
    return tiA;
}


// ************************************************************************* //
