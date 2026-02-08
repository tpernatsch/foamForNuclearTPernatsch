/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "helperFuns.H"
#include "sampledPatch.H"
#include "vtkSurfaceWriter.H"
#include "PatchTools.H"
#include "MeshedSurface.H"
#include "DynamicList.H"
#include "triFace.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * Functions * * * * * * * * * * * * * * * //

//- Instantiations of gather for scalar, vector and tensor
template<class scalar>
tmp<scalarField> helperFuns::gather(const scalar& value);

template<class vector>
tmp<vectorField> helperFuns::gather(const vector& value);

template<class tensor>
tmp<tensorField> helperFuns::gather(const tensor& value);


void helperFuns::writePatchFieldVTK
(
    const primitivePatch& patch,
    const scalarField& psi,
    const word& patchName,
    const word& psiName,
    const fileName& instance
)
{
    vtkSurfaceWriter writer(IOstream::streamFormat::ASCII);
        
    if (Pstream::parRun())
    {
        pointField points;
        faceList   faces;
        labelList  pointsMap;
        
        scalar mergeDim = 1e-8;
        
        PatchTools::gatherAndMerge
        (
            mergeDim,
            primitivePatch
            (
                SubList<face>(patch, patch.size()),
                patch.points()
            ),
            points,
            faces,
            pointsMap   
        );

        List<scalarField> tPsi(Pstream::nProcs());
        tPsi[Pstream::myProcNo()] = psi;
        Pstream::gatherList(tPsi);

         if (Pstream::master())
         {
            scalarField allPsi
            (
                ListListOps::combine<scalarField>
                (
                    tPsi,
                    accessOp<scalarField>()
                )
            );
            writer.write
            (
                instance,
                patchName,
                points,
                faces,
                psiName,
                allPsi,
                false
            );
         }
    }
    else
    {
        writer.write
        (
            instance,
            patchName,
            patch.points(),
            patch,
            psiName,
            psi,
            false
        );
    }
}


//- Face definitions for an icosahedron with points on the unit sphere
namespace icosahedron
{
    const scalar X = 0.52573;
    const scalar Z = 0.85065;
    const scalar N = 0.0;
    
    static const List<vector> points =
    {
      {-X,N,Z}, {X,N,Z}, {-X,N,-Z}, {X,N,-Z},
      {N,Z,X}, {N,Z,-X}, {N,-Z,X}, {N,-Z,-X},
      {Z,X,N}, {-Z,X, N}, {Z,-X,N}, {-Z,-X, N}
    };
     
    static const List<labelList> faces =
    {
      {0,4,1},{0,9,4},{9,5,4},{4,5,8},{4,8,1},
      {8,10,1},{8,3,10},{5,3,8},{5,2,3},{2,7,3},
      {7,10,3},{7,6,10},{7,11,6},{11,0,6},{0,1,6},
      {6,1,10},{9,0,11},{9,11,2},{9,2,5},{7,2,11}
    };
}


autoPtr<MeshedSurface<face> >
helperFuns::meshedUnitSphere
(
    label nDivisions
)
{
    // Create initial icosahedron
    DynamicList<vector> points(icosahedron::points);
    faceList faces(icosahedron::faces);
    
    for(int k = 0; k < nDivisions; k++)
    {
        label nFaces = 0;
        label nPts = points.size();
        faceList newFaces(4*faces.size());
        points.resize(nPts + 3*faces.size());
        forAll(faces, faceI)
        {
            const face& face = faces[faceI];
            vector& p0 = points[face[0]];
            vector& p1 = points[face[1]];
            vector& p2 = points[face[2]];
            
            // Add new midpoints projected onto the unit sphere
            vector p3 = 0.5*(p0 + p1);
            vector p4 = 0.5*(p1 + p2);
            vector p5 = 0.5*(p2 + p0);
            points[nPts] = p3 / mag(p3);
            points[nPts + 1] = p4 / mag(p4);
            points[nPts + 2] = p5 / mag(p5);
            
            // Create new faces
            newFaces[nFaces] = triFace(nPts, face[1], nPts + 1);
            newFaces[nFaces + 1] = triFace(face[0], nPts, nPts + 2);
            newFaces[nFaces + 2] = triFace(nPts + 2, nPts + 1, face[2]);
            newFaces[nFaces + 3] = triFace(nPts, nPts + 1, nPts + 2);
            nPts += 3;
            nFaces += 4;
        }
        faces = newFaces;
    }
    
    autoPtr<MeshedSurface<face> > surf
    (
        new MeshedSurface<face>
        (
            vectorField(points), 
            faceList(faces), 
            surfZoneList()
        )
    );
    
    surf->stitchFaces();
    
    return surf;
}

} // End namespace Foam

// ************************************************************************* //
