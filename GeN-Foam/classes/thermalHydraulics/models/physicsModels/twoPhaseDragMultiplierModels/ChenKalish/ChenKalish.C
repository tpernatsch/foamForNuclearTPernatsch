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

#include "ChenKalish.H"
#include "fluid.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace twoPhaseDragMultiplierModels
{
    defineTypeNameAndDebug(ChenKalish, 0);
    addToRunTimeSelectionTable
    (
        twoPhaseDragMultiplierModel, 
        ChenKalish, 
        twoPhaseDragMultiplierModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseDragMultiplierModels::ChenKalish::ChenKalish
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
    logX_
    (
        IOobject
        (
            IOobject::groupName("logX", objReg.name()),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless, 1e-9)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::twoPhaseDragMultiplierModels::ChenKalish::~ChenKalish()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::twoPhaseDragMultiplierModels::ChenKalish::correct()
{
    //- Compute log of the sqrt of the Lockhart-Martinelli parameter. 
    //  For numerical stability, values of X are limited between 7e-2 and 30
    //  just like the Kottowksi-Savatteri model
    logX_ =
        log
        ( 
            min
            (
                max
                (   
                    pow(mFluid_.thermo().mu()/oFluidPtr_->thermo().mu(), 0.1)*
                    pow
                    (
                        mFluid_.flowQuality()/
                        max
                        (
                            oFluidPtr_->flowQuality(),
                            dimensionedScalar("", dimless, 1e-69)
                        ), 
                        0.9
                    )*
                    sqrt(oFluidPtr_->thermo().rho()/mFluid_.thermo().rho()),
                    dimensionedScalar("", dimless, 7e-2)
                ),
                dimensionedScalar("", dimless, 30)
            )
        );

    this->setPhi2
    (
        exp
        (
            2.0*
            (
                0.0867*sqr(logX_)
            -   0.518*logX_
            +   1.59
            )
        )
    );
}


// ************************************************************************* //
