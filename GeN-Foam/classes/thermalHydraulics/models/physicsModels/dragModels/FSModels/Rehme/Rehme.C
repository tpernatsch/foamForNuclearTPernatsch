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

#include "Rehme.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(Rehme, 0);
    addToRunTimeSelectionTable(dragModel, Rehme, FSDragModels);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::Rehme::Rehme
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
    )
{
    scalar Np(this->get<label>("numberOfPins"));
    scalar Dp(this->get<scalar>("pinDiameter"));
    scalar Dw(this->get<scalar>("wireDiameter"));
    scalar Lw(this->get<scalar>("wireLeadLen"));
    scalar wetWrapPer(this->get<scalar>("wetWrapPerimeter"));
    
    scalar Pt(Dp+1.0444*Dw);
    scalar wetPinPer(Np*constant::mathematical::pi*(Dp+Dw));
    scalar B(sqrt(Pt/Dp) + pow((7.6*(Dp+Dw)*sqr(Pt/Dp)/Lw), 2.16));
    
    A_ = wetPinPer/(wetPinPer+wetWrapPer);
    B1_ = 64*sqrt(B);
    B2_ = 0.0816*pow(B, 0.9335);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::dragModels::Rehme::fd(const scalar& Re) const
{
    return (A_*(B1_/Re + B2_/pow(Re, 0.133)));
}

void Foam::dragModels::Rehme::correctKd(volTensorField& Kd) const
{   
    #include "calcKdFromFd.H"
}


// ************************************************************************* //
