/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2017 OpenCFD Ltd.
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

Application
    mergeOrSplitBaffles

Group
    grpMeshManipulationUtilities

Description
    Detects boundary faces that share points (baffles). Either merges them or
    duplicate the points.

Usage
    \b mergeOrSplitBaffles [OPTION]

    Options:
      - \par -detectOnly
        Detect baffles and write to faceSet duplicateFaces.

      - \par -split
        Detect baffles and duplicate the points (used so the two sides
        can move independently)

      - \par -dict \<dictionary\>
        Specify a dictionary to read actions from.


Note
    - can only handle pairwise boundary faces. So three faces using
      the same points is not handled (is illegal mesh anyway)

    - surfaces consisting of duplicate faces can be topologically split
    if the points on the interior of the surface cannot walk to all the
    cells that use them in one go.

    - Parallel operation (where duplicate face is perpendicular to a coupled
    boundary) is supported but not really tested.
    (Note that coupled faces themselves are not seen as duplicate faces)

\*---------------------------------------------------------------------------*/
#include "argList.H"
#include "Time.H"
#include "syncTools.H"
#include "faceSet.H"
#include "pointSet.H"
#include "meshTools.H"
#include "polyTopoChange.H"
#include "polyRemoveFace.H"
#include "polyModifyFace.H"
#include "indirectPrimitivePatch.H"
#include "processorPolyPatch.H"
#include "localPointRegion.H"
#include "duplicatePoints.H"
#include "ReadFields.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "processorMeshes.H"


using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void insertDuplicateMerge
(
    const polyMesh& mesh,
    const labelList& boundaryFaces,
    const labelList& duplicates,
    polyTopoChange& meshMod
);


label patchSize(const polyMesh& mesh, const labelList& patchIDs);


labelList patchFaces(const polyMesh& mesh, const labelList& patchIDs);


labelList findBaffles(const polyMesh& mesh, const labelList& boundaryFaces);


int removeBaffles(fvMesh& mesh, Time& runTime);//(int argc, char *argv[])


// ************************************************************************* //
