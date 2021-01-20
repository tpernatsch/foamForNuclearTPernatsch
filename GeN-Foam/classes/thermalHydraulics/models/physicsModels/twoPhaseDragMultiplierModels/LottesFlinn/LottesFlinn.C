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

#include "LottesFlinn.H"
#include "fluid.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace twoPhaseDragMultiplierModels
{
    defineTypeNameAndDebug(LottesFlinn, 0);
    addToRunTimeSelectionTable
    (
        twoPhaseDragMultiplierModel, 
        LottesFlinn, 
        twoPhaseDragMultiplierModels
    );
}
}

const Foam::Enum
<
    Foam::twoPhaseDragMultiplierModels::LottesFlinn::mode
>
Foam::twoPhaseDragMultiplierModels::LottesFlinn::modeNames_
(
    {
        { 
            mode::alpha, 
            "alpha" 
        },
        { 
            mode::XLM, 
            "Nguyen" 
        }
    }
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseDragMultiplierModels::LottesFlinn::LottesFlinn
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
    exp_
    (
        dict.lookupOrDefault<scalar>("exp", 2.0)
    ),
    mode_
    (
        modeNames_.get
        (
            this->lookupOrDefault<word>
            (
                "mode", 
                "alpha"
            )
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::twoPhaseDragMultiplierModels::LottesFlinn::~LottesFlinn()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::twoPhaseDragMultiplierModels::LottesFlinn::correct()
{
    volScalarField phi2
    (
        IOobject
        (
            "",
            mesh_.time().timeName(),
            this->db()
        ),
        mesh_,
        dimensionedScalar("", dimless, 1)
    );

    switch(mode_)
    {
        case mode::alpha :
        {
            forAll(phi2, i)
            {
                phi2[i] = pow
                (
                    max
                    (
                        mFluid_.normalized()[i], 
                        mFluid_.residualAlpha().value()
                    ),
                    -exp_
                );
            }
        }
        break;

        case mode::XLM :
        {
            const volScalarField& X(mFluid_.XLM());
            forAll(X, i)
            {
                phi2[i] = 
                    pow(1.0-pow(1.0+pow(X[i], 0.8), -0.378), -exp_);
            }
        }
        break;
    }

    this->setPhi2
    (
        phi2
    );
}


// ************************************************************************* //
