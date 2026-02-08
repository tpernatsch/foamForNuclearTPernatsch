/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "leastSquaresGradf.H"
    
#include "fvMesh.H"

#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"

#include "wedgeFvPatch.H"
#include "fvc.H"
#include "zeroGradientFvPatchFields.H"

#include "leastSquaresGradfVectors.H"
#include "gaussGrad.H"
#include "volMesh.H"
#include "surfaceMesh.H"
#include "GeometricField.H"
#include "extrapolatedCalculatedFvPatchField.H"
// #include "implicitContactFvPatchVectorField.H"
#include "implicitGapContactFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
 Foam::tmp
 <
     Foam::GeometricField
     <
         typename Foam::outerProduct<Foam::vector, Type>::type,
         Foam::fvPatchField,
         Foam::volMesh
     >
 >
leastSquaresGradf
 (
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf,
	  const GeometricField<Type, fvPatchField, volMesh>& vsf,
    const word& name
 ) 
 {
     typedef typename outerProduct<vector, Type>::type GradType;
 
     const fvMesh& mesh = ssf.mesh();
 
     tmp<GeometricField<GradType, fvPatchField, volMesh>> tlsGrad
     (
         new GeometricField<GradType, fvPatchField, volMesh>
         (
             IOobject
             (
                 name,
                 ssf.instance(),
                 mesh,
                 IOobject::NO_READ,
                 IOobject::NO_WRITE
             ),
             mesh,
             dimensioned<GradType>
             (
                 "zero",
                 ssf.dimensions()/dimLength,
                 Zero
             ),
             extrapolatedCalculatedFvPatchField<GradType>::typeName
         )
     );
     GeometricField<GradType, fvPatchField, volMesh>& lsGrad = tlsGrad.ref();
 
     // Get reference to least square vectors
     const leastSquaresGradfVectors& lsv = 
     leastSquaresGradfVectors::New(mesh);
 
     const surfaceVectorField& ownLs = lsv.pVectors();
     const surfaceVectorField& neiLs = lsv.nVectors();
 
     const labelUList& own = mesh.owner();
     const labelUList& nei = mesh.neighbour();
 
     forAll(own, facei)
     {
         label ownFacei = own[facei];
         label neiFacei = nei[facei];
 
         Type ownDeltaSf = ssf[facei] - vsf[ownFacei];
         Type neiDeltaSf = vsf[neiFacei] - ssf[facei];
 
         lsGrad[ownFacei] += ownLs[facei]*ownDeltaSf;
         lsGrad[neiFacei] -= neiLs[facei]*neiDeltaSf;
     }
 
     // Boundary faces
     forAll(vsf.boundaryField(), patchi)
     {
         const fvsPatchVectorField& patchOwnLs = ownLs.boundaryField()[patchi];
 
         const labelUList& faceCells =
             vsf.boundaryField()[patchi].patch().faceCells(); 
 
         const fvsPatchField<Type>& patchVsf = ssf.boundaryField()[patchi];

         forAll(patchVsf, patchFacei)
         {
             lsGrad[faceCells[patchFacei]] +=
                  patchOwnLs[patchFacei]
                 *(patchVsf[patchFacei] - vsf[faceCells[patchFacei]]);
         }
     }
 
 
     lsGrad.correctBoundaryConditions();
     fv::gaussGrad<Type>::correctBoundaryConditions(vsf, lsGrad);

    //  //- Additional correction due to the fact that implicitContact fvPatch 
    //  //-is not coupled (but the fvPatchField is!)
    //  typename GeometricField
    //  <
    //      typename outerProduct<vector, Type>::type, 
    //      fvPatchField, 
    //      volMesh
    //  >::Boundary& gGradbf = lsGrad.boundaryFieldRef();

    //  forAll(vsf.boundaryField(), patchi)
    //  {
    //      if 
    //      (  
    //         !vsf.boundaryField()[patchi].coupled()
    //         or             
    //         isType<implicitContactFvPatchVectorField>
    //         (
    //             vsf.boundaryField()[patchi]
    //         )
    //         or
    //         isType<implicitGapContactFvPatchVectorField>
    //         (
    //             vsf.boundaryField()[patchi]
    //         )
    //      )
    //      {
    //          const vectorField n
    //          (
    //              vsf.mesh().Sf().boundaryField()[patchi]
    //            / vsf.mesh().magSf().boundaryField()[patchi]
    //          );
 
    //          gGradbf[patchi] += n *
    //          (
    //              vsf.boundaryField()[patchi].snGrad()
    //            - (n & gGradbf[patchi])
    //          );
    //      }
    //   }
 
     return tlsGrad;
 }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
