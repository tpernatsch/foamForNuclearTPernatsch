/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "misesPlasticCreep.H"
#include "addToRunTimeSelectionTable.H"

#include "fvm.H"
#include "fvc.H"
#include "calculatedFvPatchFields.H"
#include "PstreamReduceOps.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(misesPlasticCreep, 0);
    addToRunTimeSelectionTable
    (
        misesPlasticity, 
        misesPlasticCreep, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::misesPlasticCreep::misesPlasticCreep
(
    const fvMesh& mesh,
    const dictionary& lawDict
)
:
    misesPlasticity(mesh, lawDict),
    creepModel_(creepModel::New(mesh, lawDict)),
    maxAverageCreep_(GREAT),
    maxMaximumCreep_(GREAT)
{
    const dictionary& controlDict(mesh_.time().controlDict());

    bool adjustTime = false;
    if (controlDict.found("adjustableTimeStep"))
    {
        adjustTime = controlDict.lookupOrDefault<bool>(
            "adjustableTimeStep", false);
    }
    else
    {
        adjustTime = controlDict.lookupOrDefault<bool>(
            "adjustTimeStep", false);
    }

    if(adjustTime)
    {
        maxAverageCreep_ = controlDict.lookupOrDefault<scalar>(
            "maxAverageCreep", GREAT);
        maxMaximumCreep_ = controlDict.lookupOrDefault<scalar>(
            "maxMaximumCreep", GREAT);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::misesPlasticCreep::~misesPlasticCreep()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::misesPlasticCreep::correct
(
    volScalarField& sigmaHyd, 
    volSymmTensorField& sigmaDev, 
    volSymmTensorField& epsilonEl,
    const labelList& addr
)
{
    //- Calculate new epsilon Creep increment
    creepModel_->correctCreep(epsilonEl, addr);

    misesPlasticity::correct(sigmaHyd, sigmaDev, epsilonEl, addr);
}


void Foam::misesPlasticCreep::correctEpsilonEl
(
    volSymmTensorField& epsilonEl,
    const labelList& addr
)
{
    misesPlasticity::correctEpsilonEl(epsilonEl, addr);
};


void Foam::misesPlasticCreep::correctAxialStrain
(
    volSymmTensorField& epsilon,
    const labelList& addr
)
{
    creepModel_->correctAxialStrain(epsilon, addr);
}


void Foam::misesPlasticCreep::correctAdditionalStrain
(
    volSymmTensorField& additionalStrain,
    const labelList& addr
)
{
    forAll(addr, addrI)
    {
        const label cellI(addr[addrI]);

        additionalStrain[cellI] = creepModel_->additionalStrain()[cellI]
                                + epsilonP_.internalField()[cellI];
    }
}


Foam::scalar Foam::misesPlasticCreep::nextDeltaT
(
    const labelList& addr
)
{
    //- Store old value of deltaT
    scalar oldDeltaT = mesh_.time().deltaT().value();

    //- Prepare new deltaT
    scalar deltaT(0);

    //- Prepare average and maximum plastic strain increment from last time step
    scalar averageDepsilonPEq(0);
    scalar maxDepsilonPEq(0);
    scalar volTot(0);

    const scalarField& DepsilonCreepEqI(creepModel_->DepsilonCreepEq());
    const scalarField& DepsilonPEq(DEpsilonPEq_());
    const scalarField& V(mesh_.V());

    //- Find maximum and average DepsilonCreepEq
    forAll(addr, addrI)
    {
        const label cellI = addr[addrI];
        averageDepsilonPEq += (DepsilonCreepEqI[cellI] + DepsilonPEq[cellI])*V[cellI];
        volTot += V[cellI];

        maxDepsilonPEq = max(maxDepsilonPEq, (DepsilonCreepEqI[cellI] + DepsilonPEq[cellI]));
    }

    reduce(volTot, sumOp<scalar>());

    reduce(averageDepsilonPEq, sumOp<scalar>());

    reduce(maxDepsilonPEq, maxOp<scalar>());

    averageDepsilonPEq /= volTot;

    //- Calculate average and maximum creep rate from last time step
    scalar averagePlasticStrainRate(averageDepsilonPEq/oldDeltaT);
    scalar maxPlasticStrainRate(maxDepsilonPEq/oldDeltaT);

    deltaT = min
    (
        maxAverageCreep_/max(averagePlasticStrainRate, SMALL), 
        maxMaximumCreep_/max(maxPlasticStrainRate, SMALL)
    );
    
    Info<< "Maximum deltaT calculated by misesPlasticCreep model: " 
    << mesh_.time().timeToUserTime(deltaT)
    << " with a maximum total plastic strain rate of " << maxPlasticStrainRate 
    << endl;

    return deltaT;
}
// ************************************************************************* //


