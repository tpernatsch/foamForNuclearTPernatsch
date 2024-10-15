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
#include "cuttingPlane.H"


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
#ifdef OPENFOAMFOUNDATION    
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
#elif OPENFOAMESI
        heatFlux_.reset
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
#endif           
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

    Info << endl;
}

void Foam::thermalSubSolver::heatFluxForMappings
(
    const surfaceScalarField& phi
)
{
    if(thermalOptDict_.lookupOrDefault<bool>("averageFieldsForMappings", false))
    {


        word patchName(thermalOptDict_.getOrDefault<word>("patchforAverage", "fuelOuter"));



        forAll(phi.boundaryField(), patchI)
        {
            if( (!isType<processorPolyPatch>(mesh_.boundaryMesh()[patchI])) and (mesh_.boundary()[patchI].name()==patchName))
            {
                const scalarField& phip = phi.boundaryField()[patchI];
                const scalarField& Tclad = T_.boundaryField()[patchI];
                const fvPatch& patch = mesh_.boundary()[patchI];
                const labelList faceCells = mesh_.boundaryMesh()[patchI].faceCells();
            
                forAll(faceCells, faceI)
                {
                    heatFlux_->internalFieldRef()[faceCells[faceI]] = phip[faceI]/patch.magSf()[faceI];
                    TCladForTH_->internalFieldRef()[faceCells[faceI]]= Tclad[faceI];
                }

                forAll(faceCells, cellI)
                {
                    point point(0,0,mesh_.C()[faceCells[cellI]][2]);
                    vector normal(0,0,1);
                    plane plane(point, normal);

                    // Define cutting plane
                    cuttingPlane cutPlane(plane, mesh_, true);

                    // Get cells indeces
                    const labelList& cutCells(cutPlane.meshCells());

                    forAll(cutCells, cellJ)
                    {
                        heatFlux_->internalFieldRef()[cutCells[cellJ]] = heatFlux_->internalFieldRef()[faceCells[cellI]];
                        TCladForTH_->internalFieldRef()[cutCells[cellJ]] = TCladForTH_->internalFieldRef()[faceCells[cellI]];
                    }
                }
            }
        }
    }
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
    }
    
    //Internal field of T, Cp
    const scalarField& TiOld = T_.oldTime().internalField();
    const scalarField& Ti = T_.internalField();
    
    const scalarField& CpiOld = mat_.Cp().oldTime().internalField();
    const scalarField& Cpi = mat_.Cp().internalField();
    
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalSubSolver::thermalSubSolver
(
    const fvMesh& mesh,
    materials& mat,
    const dictionary& thermalOptDict
)
:     
    mesh_(mesh),
    mat_(mat),
    thermalOptDict_(thermalOptDict),
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
    TCladForTH_(),
    mappedTFluid_(),
    mappedhtc_(),
    calculateEnthalpy_(false),
    enthalpy_(),
    residual_(0.0),
    initialResidual_(0.0),
    relResidual_(0),
    mapper_(nullptr)
{      
    if(      
#ifdef OPENFOAMFOUNDATION
        !thermalOptDict_.isNull()
#elif OPENFOAMESI
        !thermalOptDict_.isNullDict()
#endif  
    )
    {
        heatFluxSummary_ = 
        thermalOptDict_.lookupOrDefault<bool>("heatFluxSummary", false);
        
        calculateEnthalpy_ = 
        thermalOptDict_.lookupOrDefault<bool>("calculateEnthalpy", false);
    }

    mat_.correctThermoMechProperties(); 

    if(thermalOptDict_.lookupOrDefault<bool>("averageFieldsForMappings", false))
    {
        heatFlux_.reset
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
        TCladForTH_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "TCladForTH",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("TCladForTH", dimTemperature, 0),
                calculatedFvPatchField<scalar>::typeName
            )
        );
        mappedTFluid_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "mappedTFluid",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("mappedTFluid", dimTemperature, SMALL),
                zeroGradientFvPatchField<scalar>::typeName
            )
        );
        mappedhtc_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "mappedhtc",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("mappedhtc", dimPower/dimTemperature/dimArea, SMALL),
                zeroGradientFvPatchField<scalar>::typeName
            )
        );
    }
}    

Foam::autoPtr<Foam::thermalSubSolver>
Foam::thermalSubSolver::New
(
    const fvMesh& mesh,
    materials& mat,
    const dictionary& solverDict
)
{
    // Initialize type for thermalSolver class
    word type;

    dictionary thermalOptDict
    (
        solverDict.subOrEmptyDict("thermalSolverOptions")
    );

    solverDict.lookup("thermalSolver") >> type;
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
            thermalOptDict
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

    gradT_ = fvc::grad(T_);
}

 
bool Foam::thermalSubSolver::converged()
{

    bool converged(false);

    const dictionary& stressControl = 
    mesh_.solutionDict().subDict("stressAnalysis");

    bool useRelRes(stressControl.lookupOrDefault("useRelResT", true));
    
    scalar convergenceTolerance
    (
        stressControl.lookupOrDefault("T", 1e-6)
    );

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
            initialResidual_ < convergenceTolerance 
            && 
            relResidual_ < relConvergenceTolerance 
        )
        or
        (
            initialResidual_ < convergenceTolerance/10
        )
        or
        (
            relResidual_ < relConvergenceTolerance/10 
        )
        ; 
    }
    else
    {
        // In this case only the residual computed by OpenFOAM is accounted for 
        // the convergence monitoring.
        converged = 
        ( 
            initialResidual_ < convergenceTolerance
        ); 
    }

    return converged;
}



// ************************************************************************* //

