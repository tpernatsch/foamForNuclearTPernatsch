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

#include "NoKazimiFF.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(NoKazimiFF, 0);
    addToRunTimeSelectionTable(dragModel, NoKazimiFF, FFDragModels);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::NoKazimiFF::NoKazimiFF
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
    A_(0.0)
{
    //- Compute A_
    scalar P(dict.get<scalar>("pinPitch"));
    scalar D(dict.get<scalar>("pinDiameter"));
    A_ = 
        (4.0*constant::mathematical::pi)/
        (D*(2.0*sqrt(3.0)*sqr(P/D)-constant::mathematical::pi));

    const fluid& fluid1(FFPair.fluid1());
    const fluid& fluid2(FFPair.fluid2());
    if 
    (
        !(fluid1.isLiquid() and fluid2.isGas()) and
        !(fluid2.isLiquid() and fluid1.isGas())
    )
    {
        FatalErrorInFunction
            << "The NoKazimi model only works for liquid-gas systems. Set the "
            << "the stateOfMatter entry in "
            << "phaseProperties." << fluid1.name() << "Properties and/or "
            << "phaseProperties." << fluid2.name() << "Properties) to "
            << "distinguish between gas and liquid"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dragModels::NoKazimiFF::correctKd(volTensorField& Kd) const
{
    //const volScalarField& DhContinuous(FFPair_->DhContinuous());
    const volScalarField& magUr(FFPair_->magUr());
    const volScalarField& rhov(vapour_.rho());
    const volScalarField& alpha(vapour_.normalized());
    volScalarField iA(A_*alpha*min((1.0-alpha)/0.043, 1.0));

    forAll(mesh_.cells(), i)
    {
        tensor& Kdi(Kd[i]);
        scalar value
        (
            rhov[i]*magUr[i]*
            iA[i]*0.0025*
            (
                1.0 + 234.3*pow((1.0-sqrt(alpha[i])), 1.15)
            )   
        );

        Kdi[0] = value;
        Kdi[4] = value;
        Kdi[8] = value;
    }
}


// ************************************************************************* //
