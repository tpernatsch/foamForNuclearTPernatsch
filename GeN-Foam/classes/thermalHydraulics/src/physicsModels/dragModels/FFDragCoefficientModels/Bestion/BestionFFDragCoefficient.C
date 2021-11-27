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
#include "BestionFFDragCoefficient.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace FFDragCoefficientModels
{
    defineTypeNameAndDebug(Bestion, 0);
    addToRunTimeSelectionTable
    (
        FFDragCoefficientModel, 
        Bestion, 
        FFDragCoefficientModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FFDragCoefficientModels::Bestion::Bestion
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
            << "The Bestion model only works for liquid-gas systems. Set "
            << "the stateOfMatter entry in "
            << "phaseProperties." << fluid1.name() << "Properties and/or "
            << "phaseProperties." << fluid2.name() << "Properties) to "
            << "distinguish between gas and liquid"
            << exit(FatalError);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::FFDragCoefficientModels::Bestion::value
(
    const label& celli
) const
{
    scalar avap(vapour_.normalized()[celli]);
	const scalar Cm (0.188);
    scalar Co(1.2-0.2*sqrt(vapour_.rho()[celli]/liquid_.rho()[celli])); // Global Formula used in previous TRACE
    //scalar Co(1.2);   // Ancient Formula
    //scalar Co(1.0);   // New formula for the new TRACE
    scalar Ps( 
        (Co==1.0) ? 1.0 : min(pow(((1.0-Co*avap)*vapour_.magU()[celli]/max(1.0-avap,SMALL)-Co*liquid_.magU()[celli])/max(pair_.magUr()[celli],SMALL),2),0.05)
        );
    // Test for verification //
    //const volScalarField& a1(pair_.fluid1());
    //const volScalarField& a2(pair_.fluid2());
    //Info << " HERE mag UR -- " << pair_.magUr()[celli]<< endl;
    //Info << " HERE mag UL -- " << liquid_.magU()[celli]<< endl;
    //Info << " HERE mag UV -- " << vapour_.magU()[celli]<< endl;
    //Info << " HERE avap -*-  " << avap << endl;
    //Info << " HERE PS - " << Ps << endl;
    //Info << " **** KD **** " << 2.0*avap*pow(1.0-avap,3)*(pair_.DhDispersed()[celli]/pair_.DhContinuous()[celli])*(vapour_.rho()[celli]/liquid_.rho()[celli])*Ps/pow(Cm,2) << endl;
    //Info << " **** KD 2 **** " << (a1[celli]+a2[celli])*avap*pow(1.0-avap,3)/(pair_.DhContinuous()[celli])*(vapour_.rho()[celli])*pair_.magUr()[celli]*Ps/pow(Cm,2) << endl;

    return 2.0*avap*pow(1.0-avap,3)*(pair_.DhDispersed()[celli]/pair_.DhContinuous()[celli])*(vapour_.rho()[celli]/liquid_.rho()[celli])*Ps/pow(Cm,2);
}


// ************************************************************************* //
