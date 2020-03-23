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

#include "multiphysicsControl.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiphysicsControl, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiphysicsControl::multiphysicsControl
(
    const Time& runTime,
    fvMesh& THMesh,
    fvMesh& NMesh,
    fvMesh& TMMesh
)
:
    pimpleControl
    (
        THMesh,
        "PIMPLE"
    ),
    runTime_(runTime),
    thermalHydraulicMesh_(THMesh),
    neutronicMesh_(NMesh),
    thermoMechanicMesh_(TMMesh),
    topLevelDict_
    (
        IOobject
        (
            "fvSolution",
            runTime.system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    thermalHydraulicDict_(THMesh.solutionDict()),
    neutronicDict_(NMesh.solutionDict()),
    thermoMechanicDict_(TMMesh.solutionDict()),
    tightlyCoupled_(topLevelDict_.get<bool>("tightlyCoupled")),
    timeStepResidual_(topLevelDict_.get<scalar>("timeStepResidual")),
    maxTimeStepIterations_(topLevelDict_.get<label>("maxTimeStepIterations")),
    liquidFuel_
    (
        runTime.controlDict().lookupOrDefault<bool>("liquidFuel", false)
    )
{
    read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::multiphysicsControl::read()
{
    nCorrPIMPLE_ = topLevelDict_.get<label>("nOuterCorrectors");
    solveFlow_ = 
        runTime_.controlDict().lookupOrDefault<bool>
        (
            "solveFluidMechanics", true
        );
    solveEnergy_ =
        runTime_.controlDict().lookupOrDefault<bool>
        (
            "solveEnergy", false
        );
    solveNeutronics_ =
        runTime_.controlDict().lookupOrDefault<bool>
        (
            "solveNeutronics", false
        );
    solveThermoMechanics_ =
        runTime_.controlDict().lookupOrDefault<bool>
        (
            "solveThermoMechanics", false
        );

    return true;
}

// ************************************************************************* //
