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

#include "fissionGasRelease.H"
#include "addToRunTimeSelectionTable.H"
#include "globalFieldLists.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fissionGasRelease, 0);
    defineRunTimeSelectionTable(fissionGasRelease, dictionary);

    addToRunTimeSelectionTable
    (
        fissionGasRelease, 
        fissionGasRelease, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fissionGasRelease::fissionGasRelease
(
    const fvMesh& mesh,
    const materials& mat,
    const dictionary& fissionGasReleaseDict
)
:
    mesh_(mesh),
    mat_(mat),
    fissionGasReleaseDict_(fissionGasReleaseDict),
    fgr_(mesh.nCells(), 0.0),
    fgr0_(mesh.nCells(), 0.0),
    molFGR_(0.0),
    molHe_(0.0),
    maxLocalDeltaFgr_(GREAT),
    maxTotalDeltaFgr_(GREAT),
    maxFGR_(GREAT),
    timeSteppingCriterion_("maxFGR"),
    grainRadius_(createOrLookup<scalar>(mesh, "grainRadius", dimless, 0.0,
        zeroGradientFvPatchField<scalar>::typeName)),
    oxygenMetalRatio_(createOrLookup<scalar>(mesh, "oxygenMetalRatio", dimless, 0.0,
        zeroGradientFvPatchField<scalar>::typeName)),
    intragranularGasSwelling_
    (
        IOobject
        (
            "intragranularGasSwelling",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedSymmTensor
        (
            "intragranularGasSwelling", dimless, symmTensor::zero
        ),
        zeroGradientFvPatchField<scalar>::typeName
    ),
    intergranularGasSwelling_
    (
        IOobject
        (
            "intergranularGasSwelling",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedSymmTensor
        (
            "intergranularGasSwelling", dimless, symmTensor::zero
        ),
        zeroGradientFvPatchField<scalar>::typeName
    )
{
    // Set maximum fgr allowed per timestep (used in case of adjustable time 
    // step)
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
        if(controlDict.found("maxFGR"))
        {
            // This method is deprecated (left temporarily to be consistent
            // with old cases)
            maxFGR_ = (readScalar(controlDict.lookup("maxFGR")));
            timeSteppingCriterion_ = "maxFGR";
        }
        else
        {
            maxLocalDeltaFgr_ = controlDict.lookupOrDefault<scalar>("maxLocalFgrChange", GREAT);
            maxTotalDeltaFgr_ = controlDict.lookupOrDefault<scalar>("maxTotalFgrChange", GREAT);
            timeSteppingCriterion_ = "maxFgrChange";
        }
    }

    const labelListList& matAddrList(mat_.matAddrList());

    // Set initial grain radius (only if current value is == 0)
    forAll(mat_.materialsList(), i)
    {
        const materialModel& materialI(mat_.materialsList()[i]);
        const labelList& addr(matAddrList[i]);

        const word& materialName(materialI.name());

        if(isA<fuelMaterial>(mat_.materialsList()[i]))
        {
            const fuelMaterial& fuelMatI(refCast<const fuelMaterial>(materialI));

            // Take initial grain radius
            scalar rGrain = fuelMatI.rGrain();

            // Take initial O/M
            scalar oxygenMetalRatio = fuelMatI.oxygenMetalRatio();
            
            scalar averageGrainRadius(0);
            scalar averageOxygenMetal(0);
            scalar totalVol(0);
            const scalarField& volumes(mesh_.V());

            forAll(addr, j)
            {
                const label cellI = addr[j];

                if(grainRadius_[cellI] < VSMALL)
                {
                    grainRadius_[cellI] = rGrain;
                }

                if(oxygenMetalRatio_[cellI] < VSMALL)
                {
                    oxygenMetalRatio_[cellI] = oxygenMetalRatio;
                }
                
                averageGrainRadius = grainRadius_[cellI]*volumes[cellI];
                averageOxygenMetal = oxygenMetalRatio_[cellI]*volumes[cellI];
                totalVol = volumes[cellI];
            }

            reduce(averageGrainRadius, sumOp<scalar>());
            reduce(averageOxygenMetal, sumOp<scalar>());
            reduce(totalVol, sumOp<scalar>());

            Info << "Average grain radius for material \"" << materialName
            << "\": " << averageGrainRadius/max(totalVol,SMALL) << "m" << nl 
            << endl;

            Info << "Average grain O/M for material \"" << materialName
            << "\": " << averageOxygenMetal/max(totalVol,SMALL) << "m" << nl 
            << endl;
        }
    }

    grainRadius_.correctBoundaryConditions();
    oxygenMetalRatio_.correctBoundaryConditions();

    intragranularGasSwelling_.oldTime();
    intergranularGasSwelling_.oldTime();

    intergranularGasSwelling_.correctBoundaryConditions();
    intragranularGasSwelling_.correctBoundaryConditions();
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fissionGasRelease>
Foam::fissionGasRelease::New
(
    const fvMesh& mesh,
    const materials& mat,
    const dictionary& solverDict
)
{
    // Default pick type
    word type("none");

    // Prepare options dict
    dictionary fissionGasReleaseDict
    (
        solverDict.subOrEmptyDict("fissionGasReleaseOptions")
    );

    if (solverDict.found("fissionGasRelease"))
    {
        if (solverDict.isDict("fissionGasRelease"))
        {
            // Case 1: fissionGasRelease is a dictionary
            fissionGasReleaseDict = solverDict.subDict("fissionGasRelease");
            fissionGasReleaseDict.lookup("type") >> type;
        }
        else
        {
            // Case 2: fissionGasRelease is a word (type)
            solverDict.lookup("fissionGasRelease") >> type;
            
            if(solverDict.found("fgrOptions"))
            {
                fissionGasReleaseDict = solverDict.subDict("fgrOptions");
                    
                WarningIn("fissionGasRelease::New(const fvMesh&, materials&, const dictinoary&)")
                    << "Options dictionary 'fgrOptions' for `fissionGasRelease` class is deprecated. " 
                    << "Please use 'fissionGasReleaseOptions' instead or even better create a 'fissionGasRelease' dict where the name of the model is indicated with the keyword 'type'." << nl << endl;
            }
            // else, fissionGasReleaseDict stays as fissionGasReleaseOptions
        }
    }
    
    Info << "Selecting fissionGasRelease: " << type << endl;

    auto cstrIter
        = dictionaryConstructorTablePtr_->find(type);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("fissionGasRelease::New(const fvMesh&,\
            const materials&, const dictionary&)")
            << "Unknown fissionGasRelease type "
            << type << nl << nl
            << "Valid types are:" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    if (debug)
    {
        Info<< "Selecting fissionGasRelease solver --> "
            << type << endl;
    }

    return autoPtr<fissionGasRelease>
    (
        cstrIter()
        (
            mesh,
            mat,
            fissionGasReleaseDict
        )
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fissionGasRelease::~fissionGasRelease()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::fissionGasRelease::nextDeltaT()
{
    if(timeSteppingCriterion_ == "maxFgrChange")
    {
        // Current deltaT
        scalar oldDeltaT = mesh_.time().deltaT().value();

        // DeltaFgr field
        scalarField deltaFgr(fgr_ - fgr0_);

        // Reference to cell volumes
        const scalarField& Vi = mesh_.V();

        // Calculate total volume and weighted average total fgr
        scalar totalVol(0);
        scalar totalDeltaFgr(0);

        const labelListList& matAddrList(mat_.matAddrList());

        forAll(mat_.materialsList(), i)
        {
            if(isA<fuelMaterial>(mat_.materialsList()[i]))
            {
                const labelList& addr(matAddrList[i]);

                scalarField Vaddr(Vi, addr);
                scalarField deltaFgrAddr(deltaFgr, addr);

                totalDeltaFgr += gSum(Vaddr*deltaFgrAddr);
                totalVol += gSum(Vaddr);
            }
        }

        totalDeltaFgr = totalDeltaFgr/max(totalVol, VSMALL);

        // Current maximum local change in fgr percent
        scalar localMaxDeltaFgr = gMax(mag(fgr_ - fgr0_));

        // Scale time step based on desired change in local fgr per cent
        scalar nextDeltaTLocal(oldDeltaT*maxLocalDeltaFgr_/max(localMaxDeltaFgr, SMALL));

        // Scale time step based on desired change in total fgr per cent
        scalar nextDeltaTTotal(oldDeltaT*maxTotalDeltaFgr_/max(totalDeltaFgr, SMALL));

        Info<< "Maximum deltaT calculated by fgr model: " 
            << mesh_.time().timeToUserTime(min(nextDeltaTLocal, nextDeltaTTotal))
            << " with a current maximum local change in fgr% of " << localMaxDeltaFgr 
            << " and a current total change in fgr% of " << totalDeltaFgr 
            << endl;

        // // Possible criterion based on swelling?
        // scalarField swellingI    = tr(intergranularGasSwelling_.internalField())/3.0;
        // scalarField swellingOldI = tr(intergranularGasSwelling_.oldTime().internalField())/3.0;
        // scalarField deltaSwelling(mag(100.0*(swellingI - swellingOldI)/max(swellingOldI, SMALL)));
        // // Scale time step based on desired change in fgr per cent
        // nextDeltaT = (oldDeltaT*maxDeltaFgr_/max(gMax(deltaSwelling), SMALL));
        
        return min(nextDeltaTLocal, nextDeltaTTotal);
    }
    else
    {
        // This method is deprecated. Left temporarily only to be consistent 
        // with old cases. To be removed in the future

        // Store old value of deltaT
        scalar oldDeltaT = mesh_.time().deltaT().value();

        // Calculate the current rate of mol released
        scalar fgrRate(molFGR_/max(oldDeltaT, SMALL));

        scalar nextDeltaT(maxFGR_/max(fgrRate, SMALL));

        Info<< "Maximum deltaT calculated by fgr model: " 
            << mesh_.time().timeToUserTime(nextDeltaT)
            << " with a release rate of " << fgrRate << " mol/s"
            << endl;
    
        return nextDeltaT;
    }
}

// ************************************************************************* //
