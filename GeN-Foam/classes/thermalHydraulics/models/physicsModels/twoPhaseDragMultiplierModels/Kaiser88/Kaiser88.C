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

#include "Kaiser88.H"
#include "fluid.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace twoPhaseDragMultiplierModels
{
    defineTypeNameAndDebug(Kaiser88, 0);
    addToRunTimeSelectionTable
    (
        twoPhaseDragMultiplierModel, 
        Kaiser88, 
        twoPhaseDragMultiplierModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseDragMultiplierModels::Kaiser88::Kaiser88
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
    logSqrtX_
    (
        IOobject
        (
            IOobject::groupName("logSqrtX", objReg.name()),
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

Foam::twoPhaseDragMultiplierModels::Kaiser88::~Kaiser88()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::twoPhaseDragMultiplierModels::Kaiser88::correct()
{
    //- Compute log of the sqrt of the sqrt of the Lockhart-Martinelli 
    //  parameter. To avoid numerical problems, X is limited in the same data
    //  range as the Kottowski-Savatteri model, i.e. 7e-2 < X < 30 which
    //  translates in:
    //      sqrt(X) > sqrt(7e-2) = 0.265 
    //      sqrt(X) < sqrt(30) = 5.477
    logSqrtX_ =
        log
        ( 
            min
            (
                max
                (   
                    pow(mFluid_.thermo().mu()/oFluidPtr_->thermo().mu(), 0.05)*
                    pow
                    (
                        mFluid_.flowQuality()/
                        max
                        (
                            oFluidPtr_->flowQuality(),
                            dimensionedScalar("", dimless, 1e-69)
                        ), 
                        0.45
                    )*
                    pow
                    (
                        oFluidPtr_->thermo().rho()/mFluid_.thermo().rho(), 
                        0.25
                    ),
                    dimensionedScalar("", dimless, 0.265)
                ),
                dimensionedScalar("", dimless, 5.477)
            )
        );

    this->setPhi2
    (
        Foam::exp
        (
            2.0*
            (
                1.48
            -   1.05*logSqrtX_
            +   0.09*sqr(logSqrtX_)
            )
        )
    );
}


// ************************************************************************* //
