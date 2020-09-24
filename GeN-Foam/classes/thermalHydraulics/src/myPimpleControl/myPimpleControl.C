/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Released 2011 OpenCFD Ltd.
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Modified code Copyright (C) 2017 OpenCFD Ltd.
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

#include "myPimpleControl.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(myPimpleControl, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::myPimpleControl::read()
{
    pimpleControl::read();

    const dictionary pimpleDict(dict());

    minNCorrPIMPLE_ = 
        pimpleDict.lookupOrDefault("minOuterCorrectors", 1);
    corrPISOUntilConvergence_ = 
        pimpleDict.lookupOrDefault("correctUntilConvergence", false);

    return true;
}

bool Foam::myPimpleControl::criteriaSatisfied()
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

    const dictionary& solverDict = mesh_.solverPerformanceDict();
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

            Pair<scalar> residuals = (useFirstPISOInitialResidual) ? 
                firstPISOPrevPIMPLEResidual(solverPerfDictEntry)
            :   maxResidual(solverPerfDictEntry);

            checked = true;

            scalar relative = 0.0;
            bool relCheck = false;

            const bool absCheck =
                (residuals.last() < residualControl_[fieldi].absTol);

            if (storeIni)
            {
                residualControl_[fieldi].initialResidual = 
                    residuals.first();
            }
            else
            {
                const scalar iniRes =
                    residualControl_[fieldi].initialResidual + ROOTVSMALL;

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

            if (debug)
            {
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
    }

    return checked && achieved;
}

template<class Type>
bool Foam::myPimpleControl::firstPISOPrevPIMPLETypeResidual
(
    const entry& solverPerfDictEntry,
    Pair<scalar>& residuals
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    const word& fieldName = solverPerfDictEntry.keyword();

    if (mesh_.foundObject<fieldType>(fieldName))
    {
        const List<SolverPerformance<Type>> sp(solverPerfDictEntry.stream());

        residuals.first() = 
            mag
            (
                sp
                [
                    sp.size()
                -   nCorrPISOInPrevPIMPLE_
                -   nCorrPISOInPrevPrevPIMPLE_
                ].initialResidual()
            );
        residuals.last() = 
            mag(sp[sp.size()-nCorrPISOInPrevPIMPLE_].initialResidual());

        return true;
    }

    return false;
}

Foam::Pair<Foam::scalar> Foam::myPimpleControl::firstPISOPrevPIMPLEResidual
(
    const entry& solverPerfDictEntry
) const
{
    Pair<scalar> residuals(0.0, 0.0);

    // Check with builtin short-circuit
    const bool ok =
    (
        firstPISOPrevPIMPLETypeResidual<scalar>
        (solverPerfDictEntry, residuals)
    ||  firstPISOPrevPIMPLETypeResidual<vector>
        (solverPerfDictEntry, residuals)
    ||  firstPISOPrevPIMPLETypeResidual<sphericalTensor>
        (solverPerfDictEntry, residuals)
    ||  firstPISOPrevPIMPLETypeResidual<symmTensor>
        (solverPerfDictEntry, residuals)
    ||  firstPISOPrevPIMPLETypeResidual<tensor>
        (solverPerfDictEntry, residuals)
    );

    if (!ok && solutionControl::debug)
    {
        Info<<"No residual for " << solverPerfDictEntry.keyword()
            << " on mesh " << mesh_.name() << nl;
    }

    return residuals;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::myPimpleControl::myPimpleControl
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
    stopLoop_(false)
{
    read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- This is identical to what is found in pimpleControl. However, by re-
//  defining it here, I allow it to use the new criteriaSatisfied() functions
//  defined by myPimpleControl
bool Foam::myPimpleControl::loop()
{
    read();

    ++corr_;

    if (debug)
    {
        Info<< algorithmName_ << " loop: corr = " << corr_ << endl;
    }

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
        mesh_.data::remove("finalIteration");
        return false;
    }

    bool completed = false;
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

            mesh_.data::remove("finalIteration");
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

            mesh_.data::add("finalIteration", true);
            converged_ = true;
        }
    }
    else
    {
        if (finalIter())
        {
            mesh_.data::add("finalIteration", true);
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


// ************************************************************************* //
