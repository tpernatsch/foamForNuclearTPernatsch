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

#include "KottowskiSavatteri.H"
#include "fluid.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace twoPhaseDragMultiplierModels
{
    defineTypeNameAndDebug(KottowskiSavatteri, 0);
    addToRunTimeSelectionTable
    (
        twoPhaseDragMultiplierModel, 
        KottowskiSavatteri, 
        twoPhaseDragMultiplierModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseDragMultiplierModels::KottowskiSavatteri::KottowskiSavatteri
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
    log10X_
    (
        IOobject
        (
            IOobject::groupName("log10X", objReg.name()),
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

Foam::twoPhaseDragMultiplierModels::KottowskiSavatteri::~KottowskiSavatteri()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::twoPhaseDragMultiplierModels::KottowskiSavatteri::correct()
{
    //- Compute log of the sqrt of the Lockhart-Martinelli parameter. Recall
    //  the correlation is valid only for 7e-2 < X < 30
    log10X_ =
        log10
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
        pow
        (
            10,
            2.0*
            (
                0.1046*sqr(log10X_)
            -   0.5098*log10X_
            +   0.6252
            )
        )
    );
}


// ************************************************************************* //
