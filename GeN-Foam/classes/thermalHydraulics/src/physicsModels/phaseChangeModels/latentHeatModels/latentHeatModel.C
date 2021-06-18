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

#include "FFPair.H"
#include "latentHeatModel.H"
#include "phaseChangeModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(latentHeatModel, 0);
    defineRunTimeSelectionTable(latentHeatModel, latentHeatModels);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::latentHeatModel::latentHeatModel
(
    const phaseChangeModel& pcm,
    const dictionary& dict, 
    const objectRegistry& objReg
)
:
    IOdictionary
    (
        IOobject
        (
            "latentHeatModel",
            pcm.mesh().time().constant(),
            pcm.mesh()
        ),
        dict
    ),
    mesh_(pcm.mesh()),
    liquid_(pcm.liquid()),
    vapour_(pcm.vapour()),
    p_(mesh_.lookupObject<volScalarField>("p")),
    dmdt_(pcm.dmdt()),
    iT_(pcm.pair().iT()),
    adjust_(this->lookupOrDefault<bool>("adjust", false)),
    LSign_((pcm.pair().fluid1().isLiquid()) ? 1.0 : -1.0)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Default
Foam::scalar Foam::latentHeatModel::value(const label& celli) const
{
    return 1.0;
}

void Foam::latentHeatModel::correctField(volScalarField& L) const
{
    forAll(mesh_.cells(), i)
    {
        L[i] = LSign_*this->value(i);
    }
    L.correctBoundaryConditions();
    this->adjust(L);
}

void Foam::latentHeatModel::adjust(volScalarField& L) const
{
    if (adjust_)
    {
        if (LSign_ > 0)
        {
            L += 
                neg(dmdt_)*
                (
                    vapour_.thermo().he() 
                -   vapour_.thermo().he(p_, iT_)
                )
            -   pos(dmdt_)*
                (
                    liquid_.thermo().he() 
                -   liquid_.thermo().he(p_, iT_)
                );
        }
        else
        {
            L -= 
                pos(dmdt_)*
                (
                    vapour_.thermo().he() 
                -   vapour_.thermo().he(p_, iT_)
                )
            -   neg(dmdt_)*
                (
                    liquid_.thermo().he() 
                -   liquid_.thermo().he(p_, iT_)
                );
        }

        Info<< "L (avg min max) ="
        << " " << L.weightedAverage(mesh_.V()).value()
        << " " << min(L).value()
        << " " << max(L).value()
        << " J/kg" << endl;
    }
}

// ************************************************************************* //
