
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

#include "hydrogenTransportSolver.H"
#include "hydrogenTransportFvPatchScalarField.H"
#include "zeroGradientFvPatchFields.H"
#include "calculatedFvPatchFields.H"
#include "HashSet.H"
// #include "HSolverOneTSS.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(hydrogenTransportSolver, 0);

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::hydrogenTransportSolver::initConcentrations()
{
//     bool hasOneTss = isType<HSolverOneTSS>(*this);

    scalarField& Css = Css_.field();
    scalarField& Cpp = Cpp_.field();
    const scalarField& Ctot = Ctot_.field();
    const scalarField& Tssp = TSSp_.field();
//     const scalarField& Tssd = TSSd_.field();

    Css = 0.0;
    Cpp = 0.0;

    const labelListList& matAddrList(mat_.matAddrList());

    forAll(matAddrList, i)
    {
        if (TSSdModel_.set(i))
        {
            const labelList& addr(matAddrList[i]);
            forAll(addr, addrI)
            {
                const label cellI = addr[addrI];
                Css[cellI] = Tssp[cellI];
                Cpp[cellI] = Ctot[cellI] - Tssp[cellI];

//                 if(hasOneTss)
//                 {
//                     Css[cellI] = Tssd[cellI];
//                     Cpp[cellI] = Ctot[cellI] - Tssd[cellI];
//                 }
//                 else
//                 {
//                     Css[cellI] = Tssp[cellI];
//                     Cpp[cellI] = Ctot[cellI] - Tssp[cellI];
//                 }
            }
        }
    }

    //- Update all boundary values
    forAll(mesh_.boundaryMesh(), patchI)
    {
        fvPatchScalarField& pCss = Css_.boundaryFieldRef()[patchI];
        fvPatchScalarField& pCpp = Cpp_.boundaryFieldRef()[patchI];

        pCss = pCss.patchInternalField();
        pCpp = pCpp.patchInternalField();
    }

    Css_.correctBoundaryConditions();
    Cpp_.correctBoundaryConditions();
}


void Foam::hydrogenTransportSolver::computeCtot()
{
    Ctot_ = Cpp_ + Css_;

    // Calculate volume averages and summary
    const scalarField& V = mesh_.V().field();

    scalar sumCtot = 0;
    scalar sumV = 0;

    const labelListList& matAddrList(mat_.matAddrList());

    forAll(mat_.materialsList(), materialI)
    {
        if (TSSdModel_.set(materialI))
        {
            const labelList& addr(matAddrList[materialI]);

            forAll(addr, i)
            {
                const label cellI = addr[i];
                sumV += V[cellI];
                sumCtot += Ctot_[cellI]*V[cellI];
            }
        }
    }

    reduce(sumV, sumOp<scalar>());
    reduce(sumCtot, sumOp<scalar>());

    Info<< "Hydrogen summary: " << nl
        << tab << "total hydrogen content: " << sumCtot << nl
        << tab <<"average hydrogen content: " << sumCtot / max(sumV, SMALL)
        << " wt.ppm" << nl
        << endl;

}


Foam::scalar Foam::hydrogenTransportSolver::hydrideBoundaryFraction
(
    const scalar& T
)
{
     return -9.93e-11*pow(T,3) + 8.48e-8*pow(T, 2) - 5.73e-5*T + 0.623;
}


Foam::scalar Foam::hydrogenTransportSolver::hydrideVolFrac
(
    const scalar& T,
    const scalar& Cpp,
    const scalar& TSSd
)
{
    //- TSSd atomic fraction
    const scalar xalpha = TSSd / MH_ / ( TSSd/MH_ + (1e6 - TSSd/MZr_));

    //- boundary fraction
    const scalar xdelta = hydrideBoundaryFraction(T);

    const scalar term1 = Cpp / MH_ / ((Cpp / MH_) + ((1e6 - Cpp)/MZr_)); 

    return min(term1 / (xdelta - xalpha), 1);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hydrogenTransportSolver::hydrogenTransportSolver
(
    fvMesh& mesh,
    const materials& mat,
    const dictionary& elementTransportDict,
    const word& solverName
)
:
    regIOobject
    (
        IOobject
        (
            typeName,
            mesh.time().constant(),
            mesh
        )
    ),
    transportSolver(mesh, mat, elementTransportDict, solverName),
    Css_
    (
        IOobject
        (
            "Css",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimless, 0.0),
        hydrogenTransportFvPatchScalarField::typeName
    ),
    Cpp_
    (
        IOobject
        (
            "Cpp",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimless, 0.0),
        zeroGradientFvPatchField<scalar>::typeName
    ),
    Ctot_
    (
        IOobject
        (
            "Ctot",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        Css_ + Cpp_
    ),
    TSSpModel_(),
    TSSdModel_(),
    HDiffusionCoefficientModel_(nullptr),
    hydrideReorientationModel_(nullptr),
    T_(mesh_.lookupObject<volScalarField>("T")),
    gradT_(mesh_.lookupObject<volVectorField>("gradT")),
    diffH_
    (
        IOobject
        (
            "DH",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimArea/dimTime, 1.0),
        calculatedFvPatchField<scalar>::typeName
    ),
    soretTerm_
    (
        IOobject
        (
            "SoretTerm",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimless, 0.0),
        calculatedFvPatchField<scalar>::typeName
    ),
    TSSd_
    (
        IOobject
        (
            "TSSd",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimless, 0.0),
        calculatedFvPatchField<scalar>::typeName
    ),
    TSSp_
    (
        IOobject
        (
            "TSSp",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimless, 0.0),
        calculatedFvPatchField<scalar>::typeName
    ),
    hydrideVolFrac_
    (
        IOobject
        (
            "hydridevolumeFraction",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0.0),
        calculatedFvPatchField<scalar>::typeName
    ),
    radialHydrides_(),
    fracRHy_(),
    //- Gas constant
    Q_("Q*", dimEnergy/dimMoles, 25500),
    MH_(1.00784),
    MZr_(91.224)
{
    //- Read options for the hydrogen solver
    dictionary HDict
    (
        elementTransportDict.
        subOrEmptyDict("HSolverOptions")
    );     
    
    //- Read heat of transport
    if (HDict.found("Q*"))
    {
    #ifdef OPENFOAMFOUNDATION
        Q_.value() = HDict.lookup<scalar>("Q*");
    #elif OPENFOAMESI
        Q_.value() = HDict.get<scalar>("Q*");
    #endif
    }
        
    //- Init TSSd and TSSp for each material(cellZone) after reading subDict "materials"
    const dictionary& hydZonesDict = HDict.subDict("materials");

    TSSdModel_.setSize(mesh_.cellZones().size());
    TSSpModel_.setSize(mesh_.cellZones().size());

    label nKeys = hydZonesDict.keys().size();
    label nFoundZones = 0;
    label nFoundKeys = 0;
    wordHashSet matchedKeys(nKeys);

    forAll(mesh_.cellZones(), zoneI)
    {
        const word& zoneName = mesh_.cellZones()[zoneI].name();
        
        if (hydZonesDict.found(zoneName))
        {
            Info << tab << "Hydrogen transport enabled in cellZone " << zoneName << endl;

            const dictionary& currentZoneDict = hydZonesDict.subDict(zoneName);
            TSSdModel_.set(zoneI, TSSdModel::New(mesh_, mat, currentZoneDict));
            TSSpModel_.set(zoneI, TSSpModel::New(mesh_, mat, currentZoneDict));


            const word& keyword = hydZonesDict.lookupEntry(zoneName, false, true).keyword();

            if (!matchedKeys.found(keyword))
            {
                matchedKeys.insert(keyword);
                nFoundKeys++;
            }

            nFoundZones++;
        }
    }

    if (nFoundKeys != nKeys)
    {
        FatalIOErrorInFunction(hydZonesDict)
            << "Supplied entries don't seem to match the available cellZones."
            << endl << abort(FatalIOError);
    }

    if (nFoundZones == 0)
    {
        FatalIOErrorInFunction(hydZonesDict)
            << "No cellZones found in which hydrogen transport is enabled."
            << endl << abort(FatalIOError);
    }

    HDiffusionCoefficientModel_ = HDiffusionCoefficientModel::New(mesh, mat, HDict);

    if (mesh.foundObject<volScalarField>("sigmaCyl"))
    {
        radialHydrides_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "radialHydrides",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("", dimless, 0.0),
                zeroGradientFvPatchField<scalar>::typeName
            )
        );

        fracRHy_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "fracRHy",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("", dimless, 0.0),
                zeroGradientFvPatchField<scalar>::typeName
            )
        );

        hydrideReorientationModel_ = hydrideReorientationModel::New(mesh, HDict);
    };
}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::hydrogenTransportSolver::~hydrogenTransportSolver()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::hydrogenTransportSolver::correct()
{
    const scalarField& Cpp = Cpp_.oldTime().field();
    const scalarField& Css = Css_.oldTime().field();
    const scalarField& TSSd = TSSd_.field();
    const scalarField& TSSp = TSSp_.field();

    const labelListList& matAddrList(mat_.matAddrList());


    // update solubility thesholds
    forAll(mat_.materialsList(), i)
    {
        if (TSSdModel_.set(i))
        {
            const labelList& addr(matAddrList[i]);

            TSSdModel_[i].correctTSSd(addr);
            TSSpModel_[i].correctTSSp(addr);

            forAll(addr, i)
            {
                const label cellI = addr[i];

                if (TSSd[cellI] > TSSp[cellI])
                {
                    FatalErrorInFunction()
                        << "TSSd is greater than TSSp in cell " << cellI
                        << ". This suggests an input error. Please check."
                        << endl << abort(FatalError);
                }
            }
        }
    }

    TSSd_.correctBoundaryConditions();
    TSSp_.correctBoundaryConditions();

    forAll(mat_.materialsList(), i)
    {
        if (TSSdModel_.set(i))
        {
            const labelList& addr(matAddrList[i]);

            //- Update diffusion coefficient
            HDiffusionCoefficientModel_->correctDH(addr);

            //- Update hydride volume fraction
            forAll(addr, i)
            {
                const label cellI = addr[i];

                hydrideVolFrac_[cellI] = hydrideVolFrac(T_[cellI], Cpp[cellI], TSSd[cellI]);
            }

            //- Update hydride orientation
            if (hydrideReorientationModel_.valid())
            {
                hydrideReorientationModel_->correctReorientation(addr);
            }
        }
    }
}


Foam::scalar Foam::hydrogenTransportSolver::nextDeltaT()
{
    const dictionary& controlDict(mesh_.time().controlDict());
    
    scalar maxConcentrationChange = controlDict.lookupOrDefault<scalar>
    (
        "maxConcentrationChange", 
        GREAT
    );

    scalar nextDeltaT(0.0);
    scalar currentDeltaT = mesh_.time().deltaT().value();
    scalar deltaTMultiplier(1.0);

    dimensionedScalar maxDeltaCss = max(mag(Css_ - Css_.oldTime()));
    dimensionedScalar maxDeltaCpp = max(mag(Cpp_ - Cpp_.oldTime()));
    scalar maxDeltaC = max(maxDeltaCss, maxDeltaCpp).value();

    deltaTMultiplier = maxConcentrationChange / (maxDeltaC + SMALL);

    nextDeltaT = deltaTMultiplier * currentDeltaT;

    Info<< "Maximum deltaT calculated by relative change of hydrogen "
        << "concentration: " << mesh_.time().timeToUserTime(nextDeltaT)
        << " with a maximum change of " << maxDeltaC 
        << " wt.ppm"  << endl;
    
    return nextDeltaT;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

