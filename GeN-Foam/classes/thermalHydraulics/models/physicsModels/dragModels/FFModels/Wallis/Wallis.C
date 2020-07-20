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

#include "Wallis.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(Wallis, 0);
    addToRunTimeSelectionTable(dragModel, Wallis, FFDragModels);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::Wallis::Wallis
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
    vapour_
    (
        (FFPair.fluid1().isGas()) ? FFPair.fluid1() : FFPair.fluid2()
    ),
    liquid_
    (
        (FFPair.fluid1().isLiquid()) ? FFPair.fluid1() : FFPair.fluid2()
    ),
    DhStructure_
    (
        vapour_.mesh().lookupObject<volScalarField>("Dh")
    )
{
    const fluid& fluid1(FFPair.fluid1());
    const fluid& fluid2(FFPair.fluid2());
    if 
    (
        !(fluid1.isLiquid() and fluid2.isGas()) and
        !(fluid2.isLiquid() and fluid1.isGas())
    )
    {
        FatalErrorInFunction
            << "The Wallis model only works for liquid-gas systems. Set the "
            << "the stateOfMatter entry in "
            << "phaseProperties." << fluid1.name() << "Properties and/or "
            << "phaseProperties." << fluid2.name() << "Properties) to "
            << "distinguish between gas and liquid"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dragModels::Wallis::correctKd(volTensorField& Kd) const
{
    const volScalarField& magUr(FFPair_->magUr());
    const volScalarField& rhov(vapour_.rho());
    volScalarField alpha(vapour_/(vapour_+liquid_));

    forAll(mesh_.cells(), i)
    {
        tensor& Kdi(Kd[i]);
        scalar value
        (
            0.02*sqrt(alpha[i])/std::max(DhStructure_[i], 1e-3)*
            magUr[i]*rhov[i]*(1.0+150*(1.0-sqrt(alpha[i])))
        );

        Kdi[0] = value;
        Kdi[4] = value;
        Kdi[8] = value;
    }
}


// ************************************************************************* //
