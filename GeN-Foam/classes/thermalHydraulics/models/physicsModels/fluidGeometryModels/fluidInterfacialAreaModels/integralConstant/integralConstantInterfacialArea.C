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

#include "integralConstantInterfacialArea.H"
#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fluidInterfacialAreaModels
{
    defineTypeNameAndDebug(integralConstant, 0);
    addToRunTimeSelectionTable
    (
        fluidInterfacialAreaModel, 
        integralConstant, 
        fluidInterfacialAreaModels);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidInterfacialAreaModels::integralConstant::integralConstant
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
    intA0_
    (
        "value",
        dimArea,
        dict
    ),
    residualAlpha_
    (
        dimensionedScalar::lookupOrDefault
        (
            "residualAlpha",
            *this,
            dimless,
            1e-3
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidInterfacialAreaModels::integralConstant::~integralConstant()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> 
Foam::fluidInterfacialAreaModels::integralConstant::iA() const
{
    volScalarField twoPhaseRegion
    (
        pos0(dispersed_-residualAlpha_)*pos0(continuous_-residualAlpha_)
    );
    twoPhaseRegion = min(twoPhaseRegion, dimensionedScalar("", dimless, 1e-9));
    
    tmp<volScalarField> tiA
    (
        new volScalarField
        (
            twoPhaseRegion*(intA0_/fvc::domainIntegrate(twoPhaseRegion))
        )
    );
    
    return tiA;
}


// ************************************************************************* //
