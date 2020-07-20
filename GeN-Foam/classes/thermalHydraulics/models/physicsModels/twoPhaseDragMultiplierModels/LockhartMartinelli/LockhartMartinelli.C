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

#include "LockhartMartinelli.H"
#include "fluid.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace twoPhaseDragMultiplierModels
{
    defineTypeNameAndDebug(LockhartMartinelli, 0);
    addToRunTimeSelectionTable
    (
        twoPhaseDragMultiplierModel, 
        LockhartMartinelli, 
        twoPhaseDragMultiplierModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseDragMultiplierModels::LockhartMartinelli::LockhartMartinelli
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
    ),
    C_
    (
        dict.lookupOrDefault<scalar>("C", 20)
    ),
    X2_
    (
        IOobject
        (
            IOobject::groupName("LockhartMartinelliParameter", objReg.name()),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless, 1e-9)
    ),
    exp_
    (
        dict.lookupOrDefault<scalar>("exp", 0.2)
    )
{
    //- Just in case the user provides something remarkably stupid
    exp_ = min(exp_, 1.9);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::twoPhaseDragMultiplierModels::LockhartMartinelli::~LockhartMartinelli()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::twoPhaseDragMultiplierModels::LockhartMartinelli::correct()
{
    //- Compute Lockhart-Martinelli parameter
    X2_ = 
        max
        (
            pow(mFluid_.thermo().mu()/oFluidPtr_->thermo().mu(), exp_)*
            pow
            (
                mFluid_.flowQuality()/
                max
                (
                    oFluidPtr_->flowQuality(),
                    dimensionedScalar("", dimless, 1e-69)
                ), 
                2.0-exp_
            )*
            (oFluidPtr_->thermo().rho()/mFluid_.thermo().rho()),
            dimensionedScalar("", dimless, 1e-3)
        );

    phi2_ = 1.0 + C_/sqrt(X2_) + 1.0/X2_;

    this->limitPhi2();
}


// ************************************************************************* //
