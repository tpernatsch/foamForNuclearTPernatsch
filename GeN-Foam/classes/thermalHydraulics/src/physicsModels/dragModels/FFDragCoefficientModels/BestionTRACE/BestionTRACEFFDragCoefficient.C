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

#include "FFPair.H"
#include "BestionTRACEFFDragCoefficient.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace FFDragCoefficientModels
{
    defineTypeNameAndDebug(BestionTRACE, 0);
    addToRunTimeSelectionTable
    (
        FFDragCoefficientModel, 
        BestionTRACE, 
        FFDragCoefficientModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FFDragCoefficientModels::BestionTRACE::BestionTRACE
(
    const FFPair& pair,
    const dictionary& dict,
    const objectRegistry& objReg
)
:
    FFDragCoefficientModel
    (
        pair,
        dict,
        objReg
    ),
    vapour_
    (
        (pair.fluid1().isGas()) ? pair.fluid1() : pair.fluid2()
    ),
    liquid_
    (
        (!pair.fluid1().isGas()) ? pair.fluid1() : pair.fluid2()
    )
{
    const fluid& fluid1(pair.fluid1());
    const fluid& fluid2(pair.fluid2());
    if 
    (
        !(fluid1.isLiquid() and fluid2.isGas()) and
        !(fluid2.isLiquid() and fluid1.isGas())
    )
    {
        FatalErrorInFunction
            << "The BestionTRACE model only works for liquid-gas systems. Set "
            << "the stateOfMatter entry in "
            << "phaseProperties." << fluid1.name() << "Properties and/or "
            << "phaseProperties." << fluid2.name() << "Properties) to "
            << "distinguish between gas and liquid"
            << exit(FatalError);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::FFDragCoefficientModels::BestionTRACE::value
(
    const label& celli
) const
{
    scalar avap(vapour_.normalized()[celli]);
	const scalar Cm (0.188);
    //return avap*pow(1.0-avap,3)*(vapour_.rho()[celli])/(pair_.DhContinuous()[celli]*pow(Cm,2))*2/liquid_.rho()[celli];
    return 2*avap*pow(1.0-avap,3)*(vapour_.rho()[celli])/pow(Cm,2)/(pair_.DhContinuous()[celli])/pair_.rhoContinuous()[celli]*(pair_.DhDispersed()[celli]);
}


// ************************************************************************* //
