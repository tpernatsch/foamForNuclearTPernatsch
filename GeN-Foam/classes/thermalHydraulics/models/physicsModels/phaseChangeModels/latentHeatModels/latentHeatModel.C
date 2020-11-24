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

#include "latentHeatModel.H"
#include "fluid.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(latentHeatModel, 0);
    defineRunTimeSelectionTable(latentHeatModel, latentHeatModels);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::latentHeatModel::latentHeatModel
(
    const dictionary& dict, 
    const fluid& fluid1,
    const fluid& fluid2,
    const volScalarField& p,
    const volScalarField& dmdt,
    const volScalarField& iT
)
:
    IOdictionary
    (
        IOobject
        (
            "latentHeatModel",
            fluid1.mesh().time().constant(),
            fluid1.mesh()
        ),
        dict
    ),
    mesh_(fluid1.mesh()),
    liquid_
    (
        (fluid1.isLiquid()) ? fluid1 : fluid2
    ),
    vapour_
    (
        (fluid1.isGas()) ? fluid1 : fluid2
    ),
    p_(p),
    dmdt_(dmdt),
    iT_(iT),
    adjust_
    (
        this->lookupOrDefault<bool>("adjust", false)
    ),
    LSign_
    (
        (fluid1.isLiquid()) ? 1.0 : -1.0
    ),
    L_
    (
        IOobject
        (
            "latentHeat",
            fluid1.mesh().time().timeName(),
            fluid1.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fluid1.mesh(),
        dimensionedScalar("", dimEnergy/dimMass, 0.0),
        zeroGradientFvPatchScalarField::typeName
    )
{
    //- Yet another check on the stateOfMatter of the fluids. This is already
    //  done in phaseChangeModel.C, yet for now I am keeping it here too
    //  (even though it is 100% redundant)
    if 
    (
        (fluid1.isLiquid() and fluid2.isLiquid())
    or  (fluid1.isGas() and fluid2.isGas())
    )
    {
        FatalErrorInFunction
            << "latentHeat models required one liquid phase and one gas phase."
            << " Set the state of matter of each fluid via the stateOfMatter "
            << "keyword in phaseProperties." << fluid1.name() << " and "
            << "phaseProperties." << fluid2.name() << "." 
            << exit(FatalError);
    }   
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::latentHeatModel::adjust()
{
    if (adjust_)
    {
        if (LSign_ > 0)
        {
            L_ += 
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
            L_ -= 
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
        << " " << L_.weightedAverage(mesh_.V()).value()
        << " " << min(L_).value()
        << " " << max(L_).value()
        << " J/kg" << endl;
    }
}

// ************************************************************************* //
