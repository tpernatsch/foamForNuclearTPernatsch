/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2412                                                  |
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

#include "regionSolvers.H"
#include "solver.H"
#include "Time.H"
#include "interpolationCellPoint.H"
#include "radialBasisFunctionInterpolation.H"
#include "mergeOrSplitBaffles.H"
#include "multiPhysicsSolver.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionSolvers::regionSolvers(const Time& runTime)
:
    meshHandler_(runTime)
{
    List<Pair<word>> regionSolverNames;

    //Read mapping dict. Will be used later
    IOdictionary multiRegionCouplingDict
    (
        IOobject
        (
            "multiRegionCouplingDict",
            runTime.time().constant(),
            runTime.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    IOdictionary controlDict
    (
        IOobject
        (
            "controlDict",
            runTime.system(),
            runTime.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );


    if (controlDict.found("regionSolvers"))
    {
        const dictionary& regionSolversDict =
            controlDict.subDict("regionSolvers");

        forAllConstIter(dictionary, regionSolversDict, iter)
        {
            const word regionName(iter().keyword());
            const word solverName(iter().stream());

            regionSolverNames.append(Pair<word>(regionName, solverName));
        }
    }
    else //no backward compatibility for now
    {
        FatalIOErrorInFunction(runTime.controlDict())
                        << "regionSolvers list missing from "
                        << runTime.controlDict().name()
                        << exit(FatalIOError);
    }


    solvers_.setSize(regionSolverNames.size());
    prefixes_.setSize(regionSolverNames.size());

    string::size_type nRegionNameChars = 0;

    forAll(regionSolverNames, i)
    {
        const word& regionName = regionSolverNames[i].first();
        const word& solverName = regionSolverNames[i].second();

        Info << "Creating solver for region "<< regionName<<nl<<endl;

        if (solverName != "multiPhysicsSolver")
        {
            solvers_.set(i, solver::New(solverName, meshHandler_.returnMesh(regionName)));
        }
        else
        {
            solvers_.set(i, solver::New(solverName, meshHandler_.returnMesh("dummy")));
            Foam::solvers::multiPhysicsSolver* multiPhysicsSolverPtr = dynamic_cast<Foam::solvers::multiPhysicsSolver*>(&solvers_[i]);
            multiPhysicsSolverPtr->getMapper(meshHandler_);
            multiPhysicsSolverPtr->createSolvers(regionName);
            multiPhysicsSolverPtr = nullptr;
        }

        prefixes_[i] = regionName;
        nRegionNameChars = max(nRegionNameChars, regionName.size());
    }


    nRegionNameChars++;

    prefix0_.append(nRegionNameChars, ' ');

    forAll(regionSolverNames, i)
    {
        prefixes_[i].append(nRegionNameChars - prefixes_[i].size(), ' ');
    }

    // Create coupling fields

    meshHandler_.createCouplingFields(runTime);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionSolvers::~regionSolvers()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionSolvers::setGlobalPrefix() const
{
    Pout.prefix() = prefix0_;
}


void Foam::regionSolvers::setPrefix(const label i) const
{
    Pout.prefix() = prefixes_[i];
}


void Foam::regionSolvers::resetPrefix() const
{
    Pout.prefix() = string::null;
}


// void Foam::regionSolvers::interpolateAndMapFields(const Time& runTime)
// {

//     const IOdictionary couplingDict
//     (
//         IOobject
//         (
//             "multiRegionCouplingDict",
//             runTime.time().constant(),
//             runTime.db(),
//             IOobject::MUST_READ,
//             IOobject::NO_WRITE
//         )
//     );

//     if(couplingDict.found("interpolatedMappings"))
//     {
//         const dictionary interpolatedMappingDict(couplingDict.subDict("interpolatedMappings"));

//         forAll(regionMeshes_, regioni)
//         {
//             if (interpolatedMappingDict.found(regionMeshes_[regioni].name()))
//             {
//                 const dictionary regionToDict(interpolatedMappingDict.subDict(regionMeshes_[regioni].name()));

//                 scalarList fieldValues(0);
//                 scalarList xPos(0);
//                 scalarList yPos(0);
//                 scalarList zPos(0);

//                 List<word> regionsFrom(regionToDict.get<List<word>>("fromWhichRegions"));
//                 word fieldFromName(regionToDict.get<word>("fieldFromName"));
//                 word fieldToName(regionToDict.get<word>("fieldToName"));
//                 scalarList axialLocs(regionToDict.get<List<scalar>>("axialLocations"));
//                 word interpolationType(regionToDict.get<word>("interpolationType"));
//                 volScalarField& fieldToBeMapped(const_cast<volScalarField&>(regionMeshes_[regioni].lookupObject<volScalarField>(fieldToName)));


//                 forAll(regionsFrom, regionFromi)
//                 {
//                     label whichMesh(0);

//                     // Get mesh reference
//                     forAll(regionMeshes_, i)
//                     {
//                         if(regionMeshes_[i].name()==regionsFrom[regionFromi])
//                             whichMesh = i;
//                     }


//                     // - Get axial locations
//                     zPos.append(axialLocs);

//                     // - Get x and y
//                     vector centerOfMass(gSum(regionMeshes_[whichMesh].C().field()*regionMeshes_[whichMesh].V().field())/gSum(regionMeshes_[whichMesh].V().field()));

//                     for(label i = 0; i<axialLocs.size(); i++)
//                     {
//                         xPos.append(centerOfMass[0]);
//                         yPos.append(centerOfMass[1]);
//                         // - Now I have a list of (x,y,z) for one region. I need to get to associate a value to each coordinate

//                         point samplePoint(centerOfMass[0], centerOfMass[1], axialLocs[i]);

//                         // interpolationCellPoint<scalar> pointInterpolator(regionMeshes_[whichMesh]);

//                         label celli = regionMeshes_[whichMesh].findCell(samplePoint);

//                         const volScalarField& field(regionMeshes_[whichMesh].lookupObject<volScalarField>(fieldFromName));

//                         fieldValues.append(field[celli]);

//                     }
//                 }

//                 scalarList interpolationWeights(0);
//                 scalar interpolatedValue(0);


//                 if (interpolationType == "kriging")
//                 {
//                     Foam::radialBasisFunctionInterpolation::solveKriging
//                     (
//                         xPos, yPos, zPos, invRBFmatrix_
//                     );

//                     forAll(regionMeshes_[regioni].C(), centerI)
//                     {
//                         interpolatedValue = Foam::radialBasisFunctionInterpolation::kriging
//                         (
//                             xPos, yPos, zPos, fieldValues,
//                             regionMeshes_[regioni].C()[centerI][0], regionMeshes_[regioni].C()[centerI][1],regionMeshes_[regioni].C()[centerI][2],
//                             invRBFmatrix_
//                         );

//                         fieldToBeMapped[centerI] = interpolatedValue;
//                     }

//                     fieldToBeMapped.correctBoundaryConditions();
//                 }
//                 else if (interpolationType == "polyharmonicSpline")
//                 {
//                     interpolationWeights = Foam::radialBasisFunctionInterpolation::solvePolyharmonicSpline
//                     (
//                         xPos, yPos, zPos, fieldValues,  invRBFmatrix_
//                     );

//                     forAll(regionMeshes_[regioni].C(), centerI)
//                     {
//                         interpolatedValue = Foam::radialBasisFunctionInterpolation::polyharmonicSpline
//                         (
//                             interpolationWeights,
//                             xPos, yPos, zPos,
//                             regionMeshes_[regioni].C()[centerI][0], regionMeshes_[regioni].C()[centerI][1],regionMeshes_[regioni].C()[centerI][2]
//                         );

//                         fieldToBeMapped[centerI] = interpolatedValue;
//                     }
//                 }
//                 else
//                 {
//                     FatalErrorInFunction
//                         << interpolationType << " is an incorrect "
//                         << "radial basis function method. Available methods: polyharmonicSpline, gaussian, kriging"
//                         << exit(FatalError);
//                 }
//             }
//         }
//     }
// }

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::solver& Foam::regionSolvers::operator[](const label i)
{
    setPrefix(i);
    return solvers_[i];
}


// ************************************************************************* //
