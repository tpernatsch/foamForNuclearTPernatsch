/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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

#include "FinkLeibowitz.H"
#include "fluid.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace latentHeatModels
{
    defineTypeNameAndDebug(FinkLeibowitz, 0);
    addToRunTimeSelectionTable
    (
        latentHeatModel,
        FinkLeibowitz,
        latentHeatModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::latentHeatModels::FinkLeibowitz::
FinkLeibowitz
(
    const dictionary& dict, 
    const fluid& fluid1,
    const fluid& fluid2,
    const volScalarField& p,
    const volScalarField& dmdt,
    const volScalarField& iT
)
:
    latentHeatModel
    (
        dict,
        fluid1,
        fluid2,
        p,
        dmdt, 
        iT
    )
{
    setL();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::latentHeatModels::FinkLeibowitz::setL()
{
    forAll(L_, i)
    {
        //- Compute L based on the interfacial temperature, which is always at
        //  saturation. Limit it between the T bounds of the experimental
        //  correlation. The original ones are 371 K (i.e. the fusion
        //  temperature of solid sodium) and 2503.7 K (i.e. the critical 
        //  temperature of sodium). However, at the critical temperature the
        //  latent heat is 0 (by definition of the critical and super-critical
        //  states), so to avoid that, the correlation is restricted to a 
        //  temperature ever so slightly below that
        scalar T = min(max(iT_[i], 371), 2500);
        scalar oneMinTByTc(1.0-(T/2503.7));
        L_[i] = 
            LSign_*
            (
                393370*oneMinTByTc + 4398600*pow(oneMinTByTc, 0.29302)
            );
    }
    L_.correctBoundaryConditions();
}

void Foam::latentHeatModels::FinkLeibowitz::correct()
{
    if (mesh_.relaxField(L_.name()))
        L_.storePrevIter();

    setL();
    
    this->adjust();

    L_.relax();
}

// ************************************************************************* //
