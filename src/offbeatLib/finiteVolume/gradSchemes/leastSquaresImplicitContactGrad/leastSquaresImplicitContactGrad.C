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
 
 #include "leastSquaresImplicitContactGrad.H"
 #include "leastSquaresImplicitContactVectors.H"
 #include "gaussGrad.H"
 #include "fvMesh.H"
 #include "volMesh.H"
 #include "surfaceMesh.H"
 #include "GeometricField.H"
 #include "extrapolatedCalculatedFvPatchField.H"
 // #include "implicitContactFvPatchVectorField.H"
 #include "implicitGapContactFvPatchVectorField.H"
 
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
 Foam::fv::leastSquaresImplicitContactGrad<Type>::calcGrad
 (
     const GeometricField<Type, fvPatchField, volMesh>& vsf,
     const word& name
 ) const
 {
     typedef typename outerProduct<vector, Type>::type GradType;
 
     const fvMesh& mesh = vsf.mesh();
 
     tmp<GeometricField<GradType, fvPatchField, volMesh>> tlsGrad
     (
         new GeometricField<GradType, fvPatchField, volMesh>
         (
             IOobject
             (
                 name,
                 vsf.instance(),
                 mesh,
                 IOobject::NO_READ,
                 IOobject::NO_WRITE
             ),
             mesh,
             dimensioned<GradType>
             (
                 "zero",
                 vsf.dimensions()/dimLength,
                 Zero
             ),
             extrapolatedCalculatedFvPatchField<GradType>::typeName
         )
     );
     GeometricField<GradType, fvPatchField, volMesh>& lsGrad = tlsGrad.ref();
 
     // Get reference to least square vectors
     const leastSquaresImplicitContactVectors& lsv = 
     leastSquaresImplicitContactVectors::New(mesh);
 
     const surfaceVectorField& ownLs = lsv.pVectors();
     const surfaceVectorField& neiLs = lsv.nVectors();
 
     const labelUList& own = mesh.owner();
     const labelUList& nei = mesh.neighbour();
 
     forAll(own, facei)
     {
         label ownFacei = own[facei];
         label neiFacei = nei[facei];
 
         Type deltaVsf = vsf[neiFacei] - vsf[ownFacei];
 
         lsGrad[ownFacei] += ownLs[facei]*deltaVsf;
         lsGrad[neiFacei] -= neiLs[facei]*deltaVsf;
     }
 
     // Boundary faces
     forAll(vsf.boundaryField(), patchi)
     {
         const fvsPatchVectorField& patchOwnLs = ownLs.boundaryField()[patchi];
 
         const labelUList& faceCells =
             vsf.boundaryField()[patchi].patch().faceCells();
 
         if (vsf.boundaryField()[patchi].coupled())
         {
             const Field<Type> neiVsf
             (
                 vsf.boundaryField()[patchi].patchNeighbourField()
             );

             const fvPatchField<Type>& patchVsf = vsf.boundaryField()[patchi];
 
             forAll(neiVsf, patchFacei)
             {
                if
                (
                    isType<implicitGapContactFvPatchVectorField>
                    (
                        vsf.boundaryField()[patchi]
                    )  
                )
                {
                    //- Correction for ImplicitContact BC
                     lsGrad[faceCells[patchFacei]] +=
                          patchOwnLs[patchFacei]
                         *(patchVsf[patchFacei] - vsf[faceCells[patchFacei]]);

                }
                else
                {
                    //- standard method 
                     lsGrad[faceCells[patchFacei]] +=
                         patchOwnLs[patchFacei]
                        *(neiVsf[patchFacei] - vsf[faceCells[patchFacei]]);
                }
             }
         }
         else
         {
             const fvPatchField<Type>& patchVsf = vsf.boundaryField()[patchi];
 
             forAll(patchVsf, patchFacei)
             {
                 lsGrad[faceCells[patchFacei]] +=
                      patchOwnLs[patchFacei]
                     *(patchVsf[patchFacei] - vsf[faceCells[patchFacei]]);
             }
         }
     }
 
 
     lsGrad.correctBoundaryConditions();
     gaussGrad<Type>::correctBoundaryConditions(vsf, lsGrad);

     // //- Additional correction due to the fact that implicitContact fvPatch 
     // //-is not coupled (but the fvPatchField is!)
     // typename GeometricField
     // <
     //     typename outerProduct<vector, Type>::type, 
     //     fvPatchField, 
     //     volMesh
     // >::Boundary& gGradbf = lsGrad.boundaryFieldRef();

     // forAll(vsf.boundaryField(), patchi)
     // {
     //     if 
     //     (               
     //        isType<implicitContactFvPatchVectorField>
     //        (
     //            vsf.boundaryField()[patchi]
     //        )
     //        or
     //        isType<implicitGapContactFvPatchVectorField>
     //        (
     //            vsf.boundaryField()[patchi]
     //        )
     //     )
     //     {
     //         const vectorField n
     //         (
     //             vsf.mesh().Sf().boundaryField()[patchi]
     //           / vsf.mesh().magSf().boundaryField()[patchi]
     //         );

     //         gGradbf[patchi] += n *
     //         (
     //             vsf.boundaryField()[patchi].snGrad()
     //           - (n & gGradbf[patchi])
     //         );
     //     }
     //  }
 
     return tlsGrad;
 }
 
 
 // ************************************************************************* //