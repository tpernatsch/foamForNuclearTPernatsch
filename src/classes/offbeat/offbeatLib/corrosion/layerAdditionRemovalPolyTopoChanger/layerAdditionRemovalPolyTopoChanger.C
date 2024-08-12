/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "layerAdditionRemovalPolyTopoChanger.H"
#include "polyMesh.H"
#include "polyTopoChange.H"
#include "Time.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(layerAdditionRemovalPolyTopoChanger, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::layerAdditionRemovalPolyTopoChanger::layerAdditionRemovalPolyTopoChanger
(
    const IOobject& io,
    polyMesh& mesh
)
:
    polyTopoChanger(io, mesh)
{}


Foam::layerAdditionRemovalPolyTopoChanger::layerAdditionRemovalPolyTopoChanger
(
    polyMesh& mesh,
    const Switch& thicknessFromVolume
)
:
    polyTopoChanger(mesh),
    thicknessFromVolume_(thicknessFromVolume)
{
    //- Reset modifiers PtrList
    PtrList<polyMeshModifier>& modifiers = *this;
    modifiers.clear();

    //- Initialize index for modifiers
    label modifierIndex(0);

    //- Initialize index for faceZones
    label nFaceZones(mesh_.faceZones().size());

    //- Create one modifier per patch
    forAll(mesh_.boundaryMesh(), patchI)
    {
        const polyPatch& p(mesh_.boundaryMesh()[patchI]);

        // TODO: 
        // -check what happens in restart
        // -all patches or not?
        // -check what happens with large deformation (unwanted layer addition)?

        //- Avoid setting meshModifier for patches with no faces (e.g. 
        //  defaultFaces)
        if( p.size() )
        {
            modifiers.setSize(modifierIndex + 1);

            if( mesh_.faceZones().findZoneID(p.name()) == -1 )
            {
                //- If facezone with p.name is not in faceZone list, create new
                //  faceZone

                labelList zoneList(p.size());
                boolList flipList(p.size(), true);

                forAll(p, facei)
                {                
                    zoneList[facei] = p.start()+facei;
                }

                mesh_.faceZones().setSize(nFaceZones+1);
                mesh_.faceZones().set
                (
                    nFaceZones,
                    new faceZone
                    (
                        p.name(),
                        zoneList,
                        flipList,
                        nFaceZones,
                        mesh_.faceZones()
                    )
                );

                nFaceZones++;

                vectorField faceCentres(p.faceCentres());
                vectorField faceCellCentres(p.faceCellCentres());
                vectorField deltas(faceCellCentres - faceCentres);

                //- 1/4 face delta = 1/2 initial cell thickness
                scalar minThickness(gMin(mag(deltas))/4.0);

                //- 4 times face delta = 2 times initial cell thickness
                scalar maxThickness(gMax(mag(deltas))*4.0);

                modifiers.set
                (
                    modifierIndex,
                    new layerAdditionRemoval
                    (
                        p.name(),
                        modifierIndex,
                        *this,
                        p.name(),
                        minThickness,
                        maxThickness,
                        thicknessFromVolume_
                    )
                );

                modifierIndex++;
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //


// ************************************************************************* //