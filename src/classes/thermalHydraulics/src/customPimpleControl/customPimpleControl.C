/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2312                                                  |
|    Copyright 2011-2016 OpenFOAM Foundation, 2017-2025 OpenCFD Ltd.          |
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

#include "customPimpleControl.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(customPimpleControl, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::customPimpleControl::read()
{
    pimpleControl::read();

    const dictionary pimpleDict(dict());

    minNCorrPIMPLE_ =
        pimpleDict.lookupOrDefault("minNOuterCorrectors", 1);
    corrPISOUntilConvergence_ =
        pimpleDict.lookupOrDefault("correctUntilConvergence", false);

    return true;
}

bool Foam::customPimpleControl::criteriaSatisfied()
{
    // no checks on first iteration - nothing has been calculated yet
    if
    (
        (corr_ == 1)
    ||  residualControl_.empty()
    ||  finalIter()
    ||  corr_ < minNCorrPIMPLE_
    )
    {
        return false;
    }

    const bool storeIni = this->storeInitialResiduals();

    bool achieved = true;
    bool checked = false;    // safety that some checks were indeed performed

    const dictionary& solverDict = mesh_.data().solverPerformanceDict();
    forAllConstIters(solverDict, iter)
    {
        const entry& solverPerfDictEntry = *iter;

        const word& fieldName = solverPerfDictEntry.keyword();
        const label fieldi = applyToField(fieldName);

        if (fieldi != -1)
        {
            bool useFirstPISOInitialResidual
            (
                dict().subDict
                (
                    "residualControl"
                ).subDict
                (
                    fieldName
                ).lookupOrDefault<bool>("useFirstPISOInitialResidual", false)
            );

            bool stopIfFirstResRelChangeAboveRelTol
            (
                dict().subDict
                (
                    "residualControl"
                ).subDict
                (
                    fieldName
                ).lookupOrDefault<bool>
                (
                    "stopIfFirstResRelChangeAboveRelTol", false
                )
                and useFirstPISOInitialResidual
            );

            Pair<scalar> residuals = (useFirstPISOInitialResidual)
                ? previousPIMPLEFirstPISORes_
                : maxResidual(solverPerfDictEntry);

            checked = true;

            scalar relative = 0.0;
            bool relCheck = false;

            const bool absCheck =
                (residuals.last() < residualControl_[fieldi].absTol);

            if (storeIni)
            {
                residualControl_[fieldi].initialResidual = residuals.first();
            }
            else
            {
                const scalar iniRes =
                    residualControl_[fieldi].initialResidual + SMALL;

                relative =
                    (stopIfFirstResRelChangeAboveRelTol) ?
                    residuals.last() / residuals.first() :
                    residuals.last() / iniRes;

                relCheck =
                    (stopIfFirstResRelChangeAboveRelTol) ?
                    (relative >= residualControl_[fieldi].relTol) :
                    (relative < residualControl_[fieldi].relTol);

                if (stopIfFirstResRelChangeAboveRelTol and relCheck)
                    stopLoop_ = true;
            }

            achieved = achieved && (absCheck || relCheck || stopLoop_);

            Info<< algorithmName_ << " loop:" << endl;

            Info<< "    " << fieldName
                << " PIMPLE iter " << corr_
                << ": ini res = "
                << residualControl_[fieldi].initialResidual
                << ", abs tol = " << residuals.last()
                << " (" << residualControl_[fieldi].absTol << ")"
                << ", rel tol = " << relative
                << " (" << residualControl_[fieldi].relTol << ")"
                << endl;
        }
    }

    return checked && achieved;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::customPimpleControl::customPimpleControl
(
    fvMesh& mesh,
    const word& dictName,
    const bool verbose
)
:
    pimpleControl(mesh, dictName, verbose),
    minNCorrPIMPLE_(1),
    corrPISOUntilConvergence_(false),
    nCorrPISOInPrevPIMPLE_(0),
    nCorrPISOInPrevPrevPIMPLE_(0),
    stopLoop_(false),
    previousPIMPLEFirstPISORes_(0, 0),
    corr_(0)
{
    read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- This is identical to what is found in pimpleControl. However, by re-
//  defining it here, I allow it to use the new criteriaSatisfied() functions
//  defined by customPimpleControl
bool Foam::customPimpleControl::loop()
{
    read();

    ++corr_;

    Info<< algorithmName_ << " loop: corr = " << corr_ << endl;

    setFirstIterFlag();

    if (corr_ == nCorrPIMPLE_ + 1)
    {
        if (!residualControl_.empty() && (nCorrPIMPLE_ != 1))
        {
            Info<< algorithmName_ << ": not converged within "
                << nCorrPIMPLE_ << " iterations" << endl;
        }

        corr_ = 0;
        nCorrPISOInPrevPrevPIMPLE_ = 0;
        stopLoop_ = false;
        mesh_.data().setFinalIteration(false);
        return false;
    }

    bool completed = false;
//  bool converged_=false;
    if (converged_ || criteriaSatisfied())
    {

        if (converged_)
        {
            if (!stopLoop_)
                Info<< algorithmName_ << ": converged in " << corr_ - 1
                    << " iterations" << endl;
            else
                Info<< algorithmName_ << ": loop interrupted due to poor "
                    << "first PISO initialResidual convergence" << endl;

            mesh_.data().setFinalIteration(false);
            corr_ = 0;
            nCorrPISOInPrevPrevPIMPLE_ = 0;
            stopLoop_ = false;
            converged_ = false;

            completed = true;
        }
        else
        {
            Info<< algorithmName_ << ": iteration " << corr_ << endl;
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
            Info<< algorithmName_ << ": iteration " << corr_ << endl;
            storePrevIterFields();
            completed = false;
        }
    }

    nCorrPISOInPrevPrevPIMPLE_ = nCorrPISOInPrevPIMPLE_;

    return !completed;
}

void Foam::customPimpleControl::updateFirstPISOPrevPIMPLE()
{
    if (corrPISO_==1)
    {
        const dictionary& solverDict = mesh_.data().solverPerformanceDict();
    
        forAllConstIters(solverDict, iter)
        {
            const entry& solverPerfDictEntry = *iter;
    
            const word& fieldName = solverPerfDictEntry.keyword();
    
            if(fieldName == "p_rgh")
            {
                const List<SolverPerformance<scalar>> sp(solverPerfDictEntry.stream());
    
                previousPIMPLEFirstPISORes_[0] = previousPIMPLEFirstPISORes_[1];
    
                previousPIMPLEFirstPISORes_[1] =    
                (
                    mag
                    (
                        sp
                        [
                            sp.size() -1
                        ].initialResidual()
                    )
                );
            }

            else
            {
                Info <<"No residual for p_rgh on mesh " << mesh_.name() << nl;
            }
        }
    }
}


// ************************************************************************* //