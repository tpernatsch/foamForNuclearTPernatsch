/*---------------------------------------------------------------------------*\
   =========                 |
   \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
    \\    /   O peration     |
     \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
 
 #include "leastSquaresGradfVectors.H"
 #include "volFields.H"
#include "regionCoupledOFFBEATFvPatch.H"
#include "symmetryFvPatch.H"
#include "wedgeFvPatch.H"
 
 // * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
 
 namespace Foam
 {
     defineTypeNameAndDebug(leastSquaresGradfVectors, 0);
 }
 
 
 // * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //
 
 Foam::leastSquaresGradfVectors::leastSquaresGradfVectors(const fvMesh& mesh)
 :
     MeshObject<fvMesh, Foam::MoveableMeshObject, leastSquaresGradfVectors>(mesh),
     pVectors_
     (
         IOobject
         (
             "LeastSquaresP",
             mesh_.pointsInstance(),
             mesh_,
             IOobject::NO_READ,
             IOobject::NO_WRITE,
             false
         ),
         mesh_,
         dimensionedVector("zero", dimless/dimLength, Zero)
     ),
     nVectors_
     (
         IOobject
         (
             "LeastSquaresN",
             mesh_.pointsInstance(),
             mesh_,
             IOobject::NO_READ,
             IOobject::NO_WRITE,
             false
         ),
         mesh_,
         dimensionedVector("zero", dimless/dimLength, Zero)
     )
 {
     calcleastSquaresGradfVectors();
 }
 
 
 // * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //
 
 Foam::leastSquaresGradfVectors::~leastSquaresGradfVectors()
 {}
 
 
 // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
 
 void Foam::leastSquaresGradfVectors::calcleastSquaresGradfVectors()
 {
     if (debug)
     {
         InfoInFunction << "Calculating least square gradient vectors" << endl;
     }
 
     const fvMesh& mesh = mesh_;
 
     // Set local references to mesh data
     const labelUList& owner = mesh_.owner();
     const labelUList& neighbour = mesh_.neighbour();
 
     const volVectorField& C = mesh.C();
     // surfaceVectorField nf = mesh.Sf()/mesh.magSf();
     const surfaceScalarField& w = mesh.weights();
 
 
     // Set up temporary storage for the dd tensor (before inversion)
     symmTensorField dd(mesh_.nCells(), Zero);
 
     forAll(owner, facei)
     {
         label own = owner[facei];
         label nei = neighbour[facei];
 
         vector ownD = (C[nei] - C[own])*(1 - w[facei]);
         vector neiD = (C[nei] - C[own])*(w[facei]);
         symmTensor ownWdd = (1/magSqr(ownD))*sqr(ownD);
         symmTensor neiWdd = (1/magSqr(neiD))*sqr(neiD);
 
         dd[own] += ownWdd;
         dd[nei] += neiWdd;
     }
 
 
     surfaceVectorField::Boundary& pVectorsBf =
         pVectors_.boundaryFieldRef();
 
     forAll(pVectorsBf, patchi)
     {
         const fvsPatchScalarField& pw = w.boundaryField()[patchi];
 
         const fvPatch& p = pw.patch();
         const vectorField pCf = p.Cf();
         const labelUList& faceCells = p.patch().faceCells();
 
         // Build the d-vectors
         vectorField pd(p.delta());

         if(pw.coupled())
         {
             const scalarField& wp = mesh.weights().boundaryField()[patchi];

             forAll(pd, patchFacei)
             {
                // const vector d = pCf[patchFacei] - C[faceCells[patchFacei]];
                const vector& d = pd[patchFacei]*(1 - wp[patchFacei]);

                 dd[faceCells[patchFacei]] +=
                    (1/magSqr(d))*sqr(d);    
             }
         }
         else 
         {
             forAll(pd, patchFacei)
             {
                 const vector& d = pd[patchFacei];
            
                 dd[faceCells[patchFacei]] +=
                    (1/magSqr(d))*sqr(d);    
             }
         }
     }
 
 
     // Invert the dd tensor
     const symmTensorField invDd(inv(dd));
 
 
     // Revisit all faces and calculate the pVectors_ and nVectors_ vectors
     forAll(owner, facei)
     {
         label own = owner[facei];
         label nei = neighbour[facei];
 
         vector ownD = (C[nei] - C[own])*(1 - w[facei]);
         vector neiD = (C[nei] - C[own])*(w[facei]);
         scalar ownMagSfByMagSqrd = 1/magSqr(ownD);
         scalar neiMagSfByMagSqrd = 1/magSqr(neiD);
 
         pVectors_[facei] = ownMagSfByMagSqrd*(invDd[own] & ownD);
         nVectors_[facei] = -neiMagSfByMagSqrd*(invDd[nei] & neiD);
     }
 
     forAll(pVectorsBf, patchi)
     {
         fvsPatchVectorField& patchLsP = pVectorsBf[patchi];
 
         const fvsPatchScalarField& pw = w.boundaryField()[patchi];
 
         const fvPatch& p = pw.patch();
         const vectorField pCf = p.Cf();
         const labelUList& faceCells = p.faceCells();
 
         // Build the d-vectors
         vectorField pd(p.delta());

         if(pw.coupled())
         {
             const scalarField& wp = mesh.weights().boundaryField()[patchi];
             
             forAll(pd, patchFacei)
             {
                // const vector d = pCf[patchFacei] - C[faceCells[patchFacei]];
                const vector& d = pd[patchFacei]*(1 - wp[patchFacei]);
 
                patchLsP[patchFacei] =
                     (1/magSqr(d))
                    *(invDd[faceCells[patchFacei]] & d);    
             }
         }
         else
         {
             forAll(pd, patchFacei)
             {
                 const vector& d = pd[patchFacei];
 
                 patchLsP[patchFacei] =
                     (1.0/magSqr(d))
                    *(invDd[faceCells[patchFacei]] & d);
             }
         }

     }
 
     if (debug)
     {
         InfoInFunction
             << "Finished calculating least square gradient vectors" << endl;
     }
 }
 
 
 bool Foam::leastSquaresGradfVectors::movePoints()
 {
     calcleastSquaresGradfVectors();
     return true;
 }
 
 
 // ************************************************************************* //