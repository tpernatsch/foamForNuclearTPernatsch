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

#include "corrosionByPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(corrosionByPatch, 0);
    addToRunTimeSelectionTable
    (
        corrosion,
        corrosionByPatch,
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::corrosionByPatch::corrosionByPatch
(
    fvMesh& mesh,
    const materials& mat,
    const dictionary& corrosionDict
)
:
    corrosion(mesh, mat, corrosionDict),
    polyTopoChanger_(mesh),
    corrosionModels_(),
    patchIDList_(),
    currentTimeIndex_(-1),
    nOuterIter_(0)
{
    const dictionary& patchDict(corrosionDict.subDict("patches"));

    wordList names(patchDict.sortedToc());

    corrosionModels_.setSize(names.size());
    patchIDList_.setSize(names.size());

    forAll(names, i)
    {
        word name = names[i];

        label patchID(mesh_.boundary().findPatchID(name));

        if( patchID == -1)
        {
            FatalErrorInFunction()
            << "The name \"" << name << "\" given in the \"patches\" subDict "
            << "in the corrosionOptions dict" << nl
            << "does not correspond to any registered patch name." << nl << nl
            << "Possible patch are:" << nl << nl;
            forAll(mesh_.boundary(), j)
            {
                FatalErrorInFunction() << mesh_.boundaryMesh()[j].name() << endl;
            }

            FatalErrorInFunction() << abort(FatalError);
        }

        // Save addressing
        patchIDList_[i] = patchID;

        // Set corrosion model
        const dictionary& currentPatchDict = patchDict.subDict(name);
        word type(currentPatchDict.lookup("corrosionModel"));

        corrosionModels_.set
        (
            i,
            corrosionModel::New(type, mesh_, currentPatchDict)
        );
    }

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::corrosionByPatch::~corrosionByPatch()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::corrosionByPatch::correct()
{
    // Decide whether to update the mesh based on the patch-specific models
    bool updateMesh(false);

    // Store oxide thickness for relaxing
    oxideThickness_.storePrevIter();
    DOxideThickness_.storePrevIter();

    // Correct models and decide whether to update the mesh or not
    forAll(corrosionModels_, i)
    {
        const label& patchID(patchIDList_[i]);
        corrosionModels_[i].correct(oxideThickness_, DOxideThickness_, patchID);

        const scalarField& S = oxideThickness_.boundaryField()[patchID];
        const scalarField& magSf = mesh_.boundary()[patchID].magSf();

        Info<< "Corrosion summary on patch "
            << mesh_.boundaryMesh()[patchID].name() << ":" << nl
            << tab << "Average oxide layer thickness: " << gSum(S*magSf)/gSum(magSf)
            << nl << endl;

        if(corrosionModels_[i].updateMesh())
        {
            updateMesh = true;
        }
    }

    oxideThickness_.relax();
    DOxideThickness_.relax();

    // Update variation of metal thickness (not identical to DOxideThickness_)
    forAll(corrosionModels_, i)
    {
        const label& patchID(patchIDList_[i]);
        corrosionModels_[i].updateDMetalThickness(
            DOxideThickness_, DMetalThickness_, patchID);
    }

    if (updateMesh && (currentTimeIndex_ != mesh_.time().timeIndex()))
    {
        currentTimeIndex_ = mesh_.time().timeIndex();

        // Do mesh changes (use inflation - put new points in topoChangeMap)
        autoPtr<mapPolyMesh> topoChangeMap = polyTopoChanger_.changeMesh(true);

        if (topoChangeMap.valid())
        {
            pointField newPoints =
            topoChangeMap().preMotionPoints();
            mesh_.movePoints(newPoints);

            mesh_.topoChanging(true);
        }
        else
        {
            mesh_.topoChanging(false);
        }

        corrosion::updateMesh();

        // Update number of outer iterations
        nOuterIter_ = 0;
    }

    nOuterIter_++;
}
// ************************************************************************* //
