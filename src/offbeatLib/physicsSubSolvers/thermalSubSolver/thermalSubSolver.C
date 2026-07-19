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

#include "thermalSubSolver.H"
#include "fvm.H"
#include "fvc.H"
#include "calculatedFvPatchFields.H"
#include "zeroGradientFvPatchFields.H"
#include "PstreamReduceOps.H"
#include "Tuple2.H"
#include "addToRunTimeSelectionTable.H"
#include "regionCoupledOFFBEATFvPatch.H"
#include "processorPolyPatch.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(thermalSubSolver, 0);
    defineRunTimeSelectionTable(thermalSubSolver, dictionary);   

    addToRunTimeSelectionTable
    (
        thermalSubSolver, 
        thermalSubSolver, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::thermalSubSolver::calcHeatFlux
(
    const surfaceScalarField& phi
)
{
    if(!heatFlux_.valid())
    { 
        heatFlux_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "heatFlux",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("heatFlux", dimPower/dimTime/dimArea, 0),
                calculatedFvPatchField<scalar>::typeName
            )
        );

        heatFluxVec_.set
        (
            new volVectorField
            (
                IOobject
                (
                    "heatFluxVec",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedVector("heatFluxVec", dimPower/dimArea, vector::zero),
                calculatedFvPatchField<scalar>::typeName
            )
        );       
    }
    
    Info << "Boundary power summary:" << endl;

    forAll(phi.boundaryField(), patchI)
    {
        if( (!isType<processorPolyPatch>(mesh_.boundaryMesh()[patchI])) )
        {
            const scalarField& phip = phi.boundaryField()[patchI];
            const fvPatch& patch = mesh_.boundary()[patchI];
        
            Info << tab << patch.name() << ": " << gSum(phip) << endl;
        
            heatFlux_->boundaryFieldRef()[patchI] = phip / patch.magSf();
        }
    }

    const volScalarField& kappa = mesh_.lookupObject<volScalarField>("k");
    heatFluxVec_() = -gradT_ * kappa;

    Info << endl;
}


void Foam::thermalSubSolver::calcEnthalpy()
{
    if(!enthalpy_.valid())
    {
        enthalpy_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "enthalpy",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("enthalpy", dimArea/(dimTime*dimTime), 0),
                zeroGradientFvPatchField<scalar>::typeName
            )
        );
        enthalpy_->storeOldTimes();
    }
    
    //Internal field of T, Cp
    const scalarField& TiOld = T_.oldTime().internalField();
    const scalarField& Ti = T_.internalField();
    
    const scalarField& CpiOld = mat_.Cp().oldTime().internalField();
    const scalarField& Cpi = mat_.Cp().internalField();
    
    const scalarField& enthalpyiOld = enthalpy_().oldTime().internalField();

    // Enthalpy computation
    const labelListList& matAddrList(mat_.matAddrList());
    forAll(mat_.materialsList(), i)// iterate over materials
    {
        if(isA<fuelMaterial>(mat_.materialsList()[i]))
        {   
            const labelList& addr(matAddrList[i]);
            
            //Loop over cells of fuelMaterial 
            forAll(addr, j)
            {
                const label cellI = addr[j];
                //- dT
                double dT  = Ti[cellI] - TiOld[cellI];
                //- Mean Cp times dT = enthalpy variation. 
                double dH = (CpiOld[cellI] + Cpi[cellI]) / 2 * dT;
                //- Add to current time the enthalpy variation
                enthalpy_()[cellI]  += dH;
            }
        }   
    }

    enthalpy_->correctBoundaryConditions();
}


void Foam::thermalSubSolver::computeGradT()
{
    // This function computes gradT correcting the mistake done by OpenFOAM at 
    // the regionCoupled boundaries. Problem is: when T is interpolated at the 
    // boundaries in case of a regionCoupled boundary the mesh weight are used
    // rather than the correct weights (that account for hGap).
    // Workaround: correct the boundary values of Tsurf with the one of T_ that
    // have been correctly calculated with the coupled boundary condition, and
    // then compute the grad with Tsurf.

    surfaceScalarField Tsurf(fvc::interpolate(T_));

    forAll(mesh_.boundary(), patchID) 
    {
        if(isType<regionCoupledOFFBEATFvPatch>(mesh_.boundary()[patchID]))
        {
            forAll(mesh_.boundary()[patchID], faceID) 
            {
                Tsurf.boundaryFieldRef()[patchID][faceID] = 
                    T_.boundaryField()[patchID][faceID];
            }
        }
    }

    gradT_ = fvc::grad(Tsurf);
}


int Foam::thermalSubSolver::nCorrectors(const word fieldName)
{
    const dictionary& stressControl =
        mesh_.solutionDict().subDict("stressAnalysis");

    int nCorrectors = 1;
    if(stressControl.found("nCorrectors"))
    {
    #if OPENFOAMFOUNDATION
        if(stressControl.isDict("nCorrectors"))
        {
            if (stressControl.subDict("nCorrectors").found(fieldName))
            {
                nCorrectors = stressControl.subDict("nCorrectors").lookup<int>(fieldName);
            }
            else if (stressControl.subDict("nCorrectors").found("default"))
            {
                nCorrectors = stressControl.subDict("nCorrectors").lookup<int>("default");
            }
            else
            {
                FatalErrorInFunction()
                << "Number of inner correctors missing for field \"" << fieldName << "\""
                << " in the fvSolution/stressAnalysis/nCorrectors sub-dict." << nl
                << "Please set this value using either the \"" << fieldName << "\"" 
                << " or \"default\" keywords" << nl << exit(FatalError);
            }
        }
        else
        {
            nCorrectors = stressControl.lookup<int>("nCorrectors");
        }
    #elif OPENFOAMESI    
        if(stressControl.isDict("nCorrectors"))
        {
            if (stressControl.subDict("nCorrectors").found(fieldName))
            {
                nCorrectors = stressControl.subDict("nCorrectors").get<int>(fieldName);
            }
            else if (stressControl.subDict("nCorrectors").found("default"))
            {
                nCorrectors = stressControl.subDict("nCorrectors").get<int>("default");
            }
            else
            {
                FatalErrorInFunction()
                << "Number of inner correctors missing for field \"" << fieldName << "\""
                << " in the fvSolution/stressAnalysis/nCorrectors sub-dict." << nl
                << "Please set this value using either the \"" << fieldName << "\"" 
                << " or \"default\" keywords" << nl << exit(FatalError);
            }
        }
        else
        {
            nCorrectors = stressControl.get<int>("nCorrectors");
        }
    #endif
    }

    return nCorrectors;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalSubSolver::thermalSubSolver
(
    const fvMesh& mesh,
    materials& mat,
    const dictionary& thermalDict
)
:     
    mesh_(mesh),
    mat_(mat),
    thermalDict_(thermalDict),
    T_
    (
        IOobject
        (
            "T",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("T", dimTemperature, 300)
    ),
    gradT_
    (
        IOobject
        (
            "gradT",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::grad(T_)
        // mesh_,
        // dimensionedVector("gradT", dimTemperature/dimLength, vector::zero)
    ),
    heatFluxSummary_(false),
    heatFlux_(),
    heatFluxVec_(),
    calculateEnthalpy_(false),
    enthalpy_(),
    residual_(0.0),
    initialResidual_(0.0),
    relResidual_(0),
    absErr_(0)
{      
    if(      
#ifdef OPENFOAMFOUNDATION
        !thermalDict_.isNull()
#elif OPENFOAMESI
        !thermalDict_.isNullDict()
#endif  
    )
    {
        heatFluxSummary_ = 
        thermalDict_.lookupOrDefault<bool>("heatFluxSummary", false);
        
        calculateEnthalpy_ = 
        thermalDict_.lookupOrDefault<bool>("calculateEnthalpy", false);
    }

    mat_.correctThermoMechProperties(); 
}    

Foam::autoPtr<Foam::thermalSubSolver>
Foam::thermalSubSolver::New
(
    const fvMesh& mesh,
    materials& mat,
    const dictionary& solverDict
)
{
    // Default pick type
    word type("constant");

    // Prepare options dict
    dictionary thermalDict
    (
        solverDict.subOrEmptyDict("thermalSolverOptions")
    );

    if (solverDict.found("thermalSolver"))
    {
        if (solverDict.isDict("thermalSolver"))
        {
            // Case 1: thermalSolver is a dictionary
            thermalDict = solverDict.subDict("thermalSolver");
            thermalDict.lookup("type") >> type;
        }
        else
        {
            // Case 2: thermalSolver is a word (type)
            solverDict.lookup("thermalSolver") >> type;
            // thermalDict stays as thermalSolverOptions
        }
    }

    // Map deprecated `fromLatestTime` into new name `constant`
    if (type == "fromLatestTime")
    {
        WarningIn("thermalSubSolver::New(const fvMesh&, materials&, const dictinoary&)")
            << "Type 'fromLatestTime' for `thermalSolver` is deprecated. " 
            << "Please use 'constant' instead." << nl << endl;
        
        type = "constant";
    }

    Info << "Selecting thermalSolver model " << type << endl;
    
    auto cstrIter
        = dictionaryConstructorTablePtr_->find(type);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("thermalSubSolver::New(fvMesh& mesh, const volScalarField& Q)")
            << "Unknown thermalSubSolver type "
            << type << nl << nl
            << "Valid types are:" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    if (debug)
    {
        Info<< "Selecting thermalSubSolver type "
            << type << endl;
    }

    return autoPtr<thermalSubSolver>(cstrIter()
        (
            mesh, 
            mat, 
            thermalDict
        ));
}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::thermalSubSolver::~thermalSubSolver()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::thermalSubSolver::correct()
{
    // Correct the thermoMech properties.
    // Necessary in case of time varying temperature field
    // (e.g. from file or with swak4Foam)
    mat_.correctThermoMechProperties();

    computeGradT();
}

 
bool Foam::thermalSubSolver::converged
(
    const scalar& residual,
    const scalar& relResidual,
    const scalar& absErr
)
{
    bool converged(false);

    const dictionary& stressControl =
    mesh_.solutionDict().subDict("stressAnalysis");

    bool useRelRes(stressControl.lookupOrDefault("useRelResT", true));

    scalar convergenceTolerance
    (
        stressControl.lookupOrDefault("T", 1e-6)
    );

    scalar absErrTolerance
    (
        stressControl.lookupOrDefault("absErrT", 0.0)
    );

    if(absErr > absErrTolerance)
    {
        if ( useRelRes )
        {
            scalar relConvergenceTolerance
            (
                readScalar(stressControl.lookup("relT"))
            );

            // Consider converged if both initialResidual and relResidual are
            // below threshold or if one of the two is 10 times below threshold.
            converged =
            (
                residual < convergenceTolerance
                &&
                relResidual < relConvergenceTolerance
            )
            or
            (
                residual < convergenceTolerance/10
            )
            or
            (
                relResidual < relConvergenceTolerance/10
            )
            ;
        }
        else
        {
            // In this case only the residual computed by OpenFOAM is accounted for
            // the convergence monitoring.
            converged =
            (
                residual < convergenceTolerance
            );
        }
    }
    else
    {
        converged = true;
    }

    return converged;
}


bool Foam::thermalSubSolver::converged()
{
    return converged(initialResidual_, relResidual_, absErr_);
}


Foam::scalar Foam::thermalSubSolver::nextDeltaT() const
{
    const dictionary& controlDict(mesh_.time().controlDict()); 

    scalar maxRelativeIncrease_ = 
        (controlDict.lookupOrDefault<scalar>(
            "maxRelativeTemperatureIncrease", GREAT));

    scalar maxRelativeDecrease_ = 
        (controlDict.lookupOrDefault<scalar>(
            "maxRelativeTemperatureDecrease", GREAT));

    //- Rate of change of relative temperature
    volScalarField dTdt = fvc::ddt(T_)/T_;

    scalar maxRate = gMax(dTdt.internalField());
    scalar minRate = gMin(dTdt.internalField());

    scalar deltaT = min(
        maxRelativeIncrease_ / max(maxRate, SMALL),
        maxRelativeDecrease_ / max(-minRate, SMALL)
    );

    Info<< "Maximum deltaT calculated by thermal subsolver: " 
        << mesh_.time().timeToUserTime(deltaT)
        << " with a relative temperature change of " 
        << max(maxRate, -minRate) << " K/s" << endl;
    
    return deltaT;
}

// ************************************************************************* //

