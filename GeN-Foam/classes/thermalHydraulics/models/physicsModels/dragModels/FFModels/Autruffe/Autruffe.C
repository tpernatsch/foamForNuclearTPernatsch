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
    vapour_
    (
        (FFPair.fluid1().isGas()) ? FFPair.fluid1() : FFPair.fluid2()
    ),
    liquid_
    (
        (FFPair.fluid1().isLiquid()) ? FFPair.fluid1() : FFPair.fluid2()
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
            << "The Autruffe model only works for liquid-gas systems. Set "
            << "the stateOfMatter entry in "
            << "phaseProperties." << fluid1.name() << "Properties and/or "
            << "phaseProperties." << fluid2.name() << "Properties) to "
            << "distinguish between gas and liquid"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dragModels::Autruffe::correctKd(volTensorField& Kd) const
{
    const volScalarField& DhContinuous(FFPair_->DhContinuous());
    const volScalarField& magUr(FFPair_->magUr());
    const volScalarField& rhov(vapour_.rho());
    const volScalarField& alpha(vapour_.normalized());

    //- Quick breakdown, generally the F-F drag models are in the form
    //  Kd = 0.5*rhoVapour*magUr*iA*f
    //  with iA being the interfacial area between phases (hardcoded by the
    //  model here) and f being the actual dimensionless drag coefficient.
    //  For Autruffe:
    //  iA = 4.31/Dh
    //  f = (1-a)*(1+75*(1-a))^0.95
    //  I actually limit f to a min of 0.005 to avoid complete decoupling
    //  between the phases

    //- I don't care about mesh boundaries, cell-by-cell is faster
    forAll(mesh_.cells(), i)
    {
        tensor& Kdi(Kd[i]);
        const scalar& a(alpha[i]);
        scalar value
        (
            (a+liquid_[i])* // Re-scale by void fraction
            2.155/(DhContinuous[i])*magUr[i]*rhov[i]*
            max
            (
                pow
                (
                    (1.0-alpha[i])*
                    (1.0+75.0*(1.0-alpha[i])),
                    0.95
                ),
                0.005
            )
        );

        Kdi[0] = value;
        Kdi[4] = value;
        Kdi[8] = value;
    }
}


// ************************************************************************* //
