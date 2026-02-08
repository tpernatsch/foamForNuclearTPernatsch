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
     
Class
    Foam::leastSquaresGradfVectors
 
Description
    Vectors for the leastSquaresGradf gradient scheme
 
SourceFiles
    leastSquaresGradfVectors.C
 
 \*---------------------------------------------------------------------------*/
 
 #ifndef leastSquaresGradfVectors_H
 #define leastSquaresGradfVectors_H
 
 #include "MeshObject.H"
 #include "fvMesh.H"
 #include "surfaceFields.H"
 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
 namespace Foam
 {
 
 /*---------------------------------------------------------------------------*\
                      Class leastSquaresGradfVectors Declaration
 \*---------------------------------------------------------------------------*/
 
 class leastSquaresGradfVectors
 :
     public MeshObject<fvMesh, MoveableMeshObject, leastSquaresGradfVectors>
 {
     // Private data
 
         //- Least-squares gradient vectors
         surfaceVectorField pVectors_;
         surfaceVectorField nVectors_;
 
 
     // Private Member Functions
 
         //- Construct Least-squares gradient vectors
         void calcleastSquaresGradfVectors();
 
 
 public:
 
     // Declare name of the class and its debug switch
     TypeName("leastSquaresGradfVectors");
 
 
     // Constructors
 
         //- Construct given an fvMesh
         explicit leastSquaresGradfVectors(const fvMesh&);
 
 
     //- Destructor
     virtual ~leastSquaresGradfVectors();
 
 
     // Member functions
 
         //- Return reference to owner least square vectors
         const surfaceVectorField& pVectors() const
         {
             return pVectors_;
         }
 
         //- Return reference to neighbour least square vectors
         const surfaceVectorField& nVectors() const
         {
             return nVectors_;
         }
 
         //- Delete the least square vectors when the mesh moves
         virtual bool movePoints();
 };
 
 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
 } // End namespace Foam
 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
 #endif