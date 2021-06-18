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
#include "SchorInterfacialArea.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace interfacialAreaModels
{
    defineTypeNameAndDebug(Schor, 0);
    addToRunTimeSelectionTable
    (
        interfacialAreaModel,
        Schor,
        interfacialAreaModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfacialAreaModels::Schor::Schor
(
    const FFPair& pair,
    const dictionary& dict,
    const objectRegistry& objReg
)
:
    interfacialAreaModel
    (
        pair,
        dict,
        objReg
    ),
    vapour_
    (
        (pair.fluid1().isGas()) ?
        pair.fluid1()
        :
        pair.fluid2()
    ),
    D_(dict.get<scalar>("pinDiameter")),
    P_(dict.get<scalar>("pinPitch")),
    PD_(P_/D_),
    A_(2.0*Foam::sqrt(3.0)*Foam::sqr(PD_)),
    PI_(3.1415927),
    alpha1_(0.55),
    alpha2_(0.65),
    deltaAlpha21_(alpha2_-alpha1_),
    alpha3_(dict.getOrDefault<scalar>("cutoffAlpha", 0.957)),
    iA1_((4.0/D_)*Foam::sqrt(PI_*alpha1_/(3.0*(A_-PI_)))),
    iA2_
    (
        (4.0/D_)*Foam::sqrt(PI_)
        /
        (
            3.0*(A_-PI_)
        )*
        Foam::sqrt
        (
            (1.0-alpha2_)*A_ + PI_*alpha2_
        )
    ),
    iA3_
    (
        dict.getOrDefault<scalar>("minInterfacialAreaAtLargeAlpha", 0)
    )
{
    if 
    (
        !(pair.fluid1().isLiquid() and pair.fluid2().isGas())
    and !(pair.fluid2().isLiquid() and pair.fluid1().isGas())
    )
    {
        FatalErrorInFunction
            << "Shor model only work for gas-liquid systems! Set the "
            << "stateOfMatter keyword in the fluid properties dictionaries to"
            << "specify the stateOfMatter of each fluid"
            << exit(FatalError);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::interfacialAreaModels::Schor::value
(
    const label& celli
) const
{
    const scalar& alpha(vapour_.normalized()[celli]);
    scalar scaleFactor(vapour_[celli]/max(alpha, 1e-9));
    return 
        scaleFactor*
        (
            (alpha < alpha1_) ?
            (4.0/D_)*Foam::sqrt(PI_*alpha/(3.0*(A_-PI_)))
            :
            (
                (alpha < alpha2_) ?
                ((alpha2_-alpha)*iA1_+(alpha-alpha1_)*iA2_)/deltaAlpha21_
                :
                max
                (
                    (4.0/D_)*Foam::sqrt(PI_)
                    /
                    (
                        3.0*(A_-PI_)
                    )*
                    Foam::sqrt
                    (
                        (1.0-alpha)*A_ + PI_*alpha
                    )*
                    min
                    (
                        (1-alpha)/(1.0-alpha3_),
                        1.0
                    ),
                    iA3_/scaleFactor
                )
            )
        );
}

// ************************************************************************* //
