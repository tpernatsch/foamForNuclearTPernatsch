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

#include "Kaiser74.H"
#include "fluid.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace twoPhaseDragMultiplierModels
{
    defineTypeNameAndDebug(Kaiser74, 0);
    addToRunTimeSelectionTable
    (
        twoPhaseDragMultiplierModel, 
        Kaiser74, 
        twoPhaseDragMultiplierModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseDragMultiplierModels::Kaiser74::Kaiser74
(
    const objectRegistry& objReg,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    twoPhaseDragMultiplierModel
    (
        objReg,
        dict,
        mesh
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::twoPhaseDragMultiplierModels::Kaiser74::~Kaiser74()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::twoPhaseDragMultiplierModels::Kaiser74::correct()
{
    /*
    So 
        phi = 8.2/X^0.55    ->
        phi2 = 67.24/(X2)^0.55
    The exponents in X2 are 0.2 for the mu ratio, 1.8 for the flow quality
    ratio and 1 for the density ratio. Multiplied by 0.55, they become the
    values below
    */
    this->setPhi2
    (
        67.24/
        max
        (   
            pow(mFluid_.thermo().mu()/oFluidPtr_->thermo().mu(), 0.11)*
            pow
            (
                mFluid_.flowQuality()/
                max
                (
                    oFluidPtr_->flowQuality(),
                    dimensionedScalar("", dimless, 1e-69)
                ), 
                0.99
            )*
            pow
            (
                oFluidPtr_->thermo().rho()/mFluid_.thermo().rho(), 
                0.55
            ),
            dimensionedScalar("", dimless, 1e-9)
        )
    );
}


// ************************************************************************* //
