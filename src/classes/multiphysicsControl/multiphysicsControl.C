/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2406                                                  |
|    Copyright 2011-2016 OpenFOAM Foundation, 2017-2024 OpenCFD Ltd.          |
-------------------------------------------------------------------------------
License
    This file is part of GeN-Foam.

    GeN-Foam is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    GeN-Foam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    This offering is not approved or endorsed by the OpenFOAM Foundation nor
    OpenCFD Limited, producer and distributor of the OpenFOAM(R)software via
    www.openfoam.com, and owner of the OPENFOAM(R) and OpenCFD(R) trademarks.

    This particular snippet of code is developed according to the developer's
    knowledge and experience in OpenFOAM. The users should be aware that
    there is a chance of bugs in the code, though we've thoroughly test it.
    The source code may not be in the OpenFOAM coding style, and it might not
    be making use of inheritance of classes to full extent.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "multiphysicsControl.H"

#if defined __has_include
#  if __has_include(<commDataLayer.H>) 
#    include <commDataLayer.H>
#    define isCommDataLayerIncluded
#  endif
#endif

#ifdef isCommDataLayerIncluded
#include "commDataLayer.H"
#endif

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiphysicsControl, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiphysicsControl::multiphysicsControl
(
    Time& runTime,
    fvMesh& THMesh,
    fvMesh& NMesh,
    fvMesh& TMMesh
)
:
    customPimpleControl
    (
        THMesh,
        runTime,
        "PIMPLE"
    ),
    // runTime_(runTime),
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
    #ifdef isCommDataLayerIncluded
    commDataLayer& data = commDataLayer::New(runTime_);
    data.storeObj
    (
        1.0,
        "isEnergyNeutronicsThermMechCorrectorConverged",
        commDataLayer::causality::out
    );
    #endif

    read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Copied straight from customPimpleControl to modify infos
bool Foam::multiphysicsControl::loop()
{
    read();
 
    ++corr_;
 
    // Extract label index of FMUController type functionObject, 
    // functionObjects are created by runTime.run()
    // Add if in first iteration, first time index
    #ifdef isCommDataLayerIncluded
    label FMUControllerLabel(-1);
    forAll(runTime_.functionObjects(), labelI)
    {
        if (runTime_.functionObjects()[labelI].type() == "FMUController")
        {
            FMUControllerLabel = labelI;
            break; // FMUController is unique
        }
    }

    // FMI check implicit step // only if multiphysics loop completed
    /*if (
        FMUControllerLabel != -1 
        && !converged_
    ) // && completed)
    {
        runTime_.functionObjects()[FMUControllerLabel].execute();

        commDataLayer& data = commDataLayer::New(runTime_);

        label isNewStep = data.getObj<label>("new_step", commDataLayer::causality::in);
        
        // nCorrPIMPLE_ = corr_+2;
        if (isNewStep == 1)
        {
            // converged_ = true;
            nCorrPIMPLE_ = corr_;
        }
        else
        {
            converged_ = false;
            nCorrPIMPLE_ = corr_+2;
        }

        Info<< "Multiphysics loop newStep: " << isNewStep << " corr=" << corr_
            << endl;
    }*/
    #endif


    setFirstIterFlag();
 
    if (corr_ == nCorrPIMPLE_ + 1)
    {
        if (!residualControl_.empty() && (nCorrPIMPLE_ != 1))
        {
            Info<< "Outer iterations not converged within "
                << nCorrPIMPLE_ << " iterations" << endl;
        }
 
        corr_ = 0;

        mesh_.data().setFinalIteration(false);

        #ifdef isCommDataLayerIncluded
        // FMI check implicit step
        if (FMUControllerLabel != -1)
        {
            runTime_.functionObjects()[FMUControllerLabel].execute();

            commDataLayer& data = commDataLayer::New(runTime_);

            label isNewStep = data.getObj<label>("new_step", commDataLayer::causality::in);
            
            Info<< "Multiphysics loop newStep last iter: " << isNewStep 
                << endl;

            return(isNewStep != 1);
        }
        #endif
        
        return false;
    }
 
    bool completed = false;
    if (converged_ || customPimpleControl::criteriaSatisfied())
    {
        if (converged_)
        {
            Info<< algorithmName_ << ": converged in " << corr_ - 1
                << " iterations" << endl;
 
            mesh_.data().setFinalIteration(false);
            corr_ = 0;
            converged_ = false;
 
            completed = true;
        }
        else
        {
            //- Neutronics and thermoMechanics are solved on the last iteration
            //  only anyway, so print Outer loop iteration info only if solving
            //  any of them
            if (solveFlow_ or solveEnergy_)
                Info<< "Outer iteration " << corr_ << endl;
            storePrevIterFields();
 
            mesh_.data().setFinalIteration(true);
            converged_ = true;
        }
    }
    else
    {
        if (finalIter())
        {
            mesh_.data().setFinalIteration(true);
        }
 
        if (corr_ <= nCorrPIMPLE_)
        {
            if (solveFlow_ or solveEnergy_)
                Info<< "Outer iteration " << corr_ << endl;
            storePrevIterFields();
            completed = false;
        }
    }

    
    #ifdef isCommDataLayerIncluded
    // FMI check implicit step // only if multiphysics loop completed
    if (FMUControllerLabel != -1 && completed)
    {
        runTime_.functionObjects()[FMUControllerLabel].execute();

        commDataLayer& data = commDataLayer::New(runTime_);

        label isNewStep = data.getObj<label>("new_step", commDataLayer::causality::in);
        
        Info<< "Multiphysics loop newStep: " << isNewStep 
            << endl;

        return(isNewStep != 1);
    }
    #endif
    

    return !completed;
}


bool Foam::multiphysicsControl::loopEnergyNeutronicsThermomechanics
(
    scalar couplingResidual,
    scalar couplingIter
)
{
    if (couplingIter == 0) // Do at least one loop at the first iteration
    {
        return(true);
    }
    
    bool isConverged = true;
    if (
        couplingResidual > timeStepResidual()
        && couplingIter < maxTimeStepIterations() 
        && tightlyCoupled()
    ) {
        isConverged = false;
    }

    Info<< "Is E-N-TM converged: " << isConverged 
        << " (isEnergyNeutroThemMechFMICorrection_: " << isEnergyNeutroThemMechFMICorrection_ << ")" 
        << endl;

    // FMI loop correction
    #ifdef isCommDataLayerIncluded
    if (isEnergyNeutroThemMechFMICorrection_) // && !(couplingIter == 1 && isConverged))
    {
        commDataLayer& data = commDataLayer::New(runTime_);
        scalar& isConvergedFMI = data.getObj<scalar>
        (
            "isEnergyNeutronicsThermMechCorrectorConverged", 
            commDataLayer::causality::out
        );
        isConvergedFMI = isConverged ? 1.0 : 0.0;
        // Need to execute the functionObj FMUController to send the info back to 
        // the main script 
        label FMUControllerLabel(-1);
        forAll(runTime_.functionObjects(), labelI)
        {
            if (runTime_.functionObjects()[labelI].type() == "FMUController")
            {
                FMUControllerLabel = labelI;
                break; // FMUController is unique
            }
        }
        if (FMUControllerLabel != -1)
        {
            runTime_.functionObjects()[FMUControllerLabel].execute();
        }
    }
    #endif

    return(!isConverged);
}


bool Foam::multiphysicsControl::read()
{
    customPimpleControl::read();
    nCorrPIMPLE_ = topLevelDict_.get<label>("nOuterCorrectors");
    isEnergyNeutroThemMechFMICorrection_ = 
        topLevelDict_.lookupOrDefault<bool>
        (
            "fmiCoupledCorrector", false
        );
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
            "solveThermalMechanics", false
        );

    return true;
}

// ************************************************************************* //
