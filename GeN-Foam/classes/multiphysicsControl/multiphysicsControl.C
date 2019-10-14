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
        "OuterLoop"
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
    thermalHydraulicDict_(topLevelDict_.subDict("thermalHydraulics")),
    neutronicDict_(topLevelDict_.subDict("neutronics")),
    thermoMechanicDict_(topLevelDict_.subDict("thermoMechanics")),
    tightlyCoupled_(topLevelDict_.get<bool>("tightlyCoupled")),
    timeStepResidual_(topLevelDict_.get<scalar>("timeStepResidual")),
    maxTimeStepIterations_(topLevelDict_.get<label>("maxTimeStepIterations")),
    integralPredictor_(neutronicDict_.get<bool>("integralPredictor")),
    implicitPredictor_(neutronicDict_.get<bool>("implicitPredictor")),
    aitkenAcceleration_(neutronicDict_.get<bool>("aitkenAcceleration")),
    neutronIterationResidual_
    (
        neutronicDict_.get<scalar>("neutronIterationResidual")
    ),
    maxNeutronIterations_(neutronicDict_.get<label>("maxNeutronIterations")),
    compactNormalStress_(thermoMechanicDict_.get<bool>("compactNormalStress")),
    nThermoMechanicCorrs_(thermoMechanicDict_.get<label>("nCorrs")),
    thermoMechanicCorr_(0),
    D_(thermoMechanicDict_.get<scalar>("D"))
{
    read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::multiphysicsControl::read()
{
    nCorrPIMPLE_ = topLevelDict_.lookupOrDefault<label>("nOuterCorrectors", 1);
    nCorrPISO_ = 
        thermalHydraulicDict_.lookupOrDefault<label>("nCorrectors", 1);
    nNonOrthCorr_ = 
        thermalHydraulicDict_.lookupOrDefault<label>("nNonOrthoCorrectors", 1);
    solveFlow_ = 
        runTime_.controlDict().lookupOrDefault<bool>
        (
            "solveFluidMechanics", true
        );
    solveEnergy_ =
        runTime_.controlDict().lookupOrDefault<bool>
        (
            "solveEnergy", true
        );
    solveNeutronics_ =
        runTime_.controlDict().lookupOrDefault<bool>
        (
            "solveNeutronics", true
        );
    solveThermoMechanics_ =
        runTime_.controlDict().lookupOrDefault<bool>
        (
            "solveThermoMechanics", true
        );
}

// ************************************************************************* //
