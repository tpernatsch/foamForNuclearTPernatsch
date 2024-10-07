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

#include "gaussLinearImplicitContactGrad.H"
#include "extrapolatedCalculatedFvPatchField.H"
#ifdef OPENFOAMFOUNDATION
// #include "implicitContactFvPatchVectorField.H"
// #include "implicitGapContactFvPatchVectorField.H"
#endif

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
Foam::fv::gaussLinearImplicitContactGrad<Type>::gradf
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf,
    const word& name
)
{
    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = ssf.mesh();

    tmp<GeometricField<GradType, fvPatchField, volMesh>> tgGrad
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
                "0",
                ssf.dimensions()/dimLength,
                Zero
            ),
            extrapolatedCalculatedFvPatchField<GradType>::typeName
        )
    );
    GeometricField<GradType, fvPatchField, volMesh>& gGrad = tgGrad.ref();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const vectorField& Sf = mesh.Sf();

    Field<GradType>& igGrad = gGrad;
    const Field<Type>& issf = ssf;

    forAll(owner, facei)
    {
        GradType Sfssf = Sf[facei]*issf[facei];

        igGrad[owner[facei]] += Sfssf;
        igGrad[neighbour[facei]] -= Sfssf;
    }

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells =
            mesh.boundary()[patchi].faceCells();

        const vectorField& pSf = mesh.Sf().boundaryField()[patchi];

        const fvsPatchField<Type>& pssf = ssf.boundaryField()[patchi];

        forAll(mesh.boundary()[patchi], facei)
        {
            igGrad[pFaceCells[facei]] += pSf[facei]*pssf[facei];
        }
    }

    igGrad /= mesh.V();

    gGrad.correctBoundaryConditions();

    return tgGrad;
}


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
Foam::fv::gaussLinearImplicitContactGrad<Type>::calcGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vsf,
    const word& name
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;

    // Create a local copy of the interpolated field, so that the values at the
    // implicit contact boundary can be modified
    GeometricField<Type, fvsPatchField, surfaceMesh> ssf = 
    tinterpScheme_().interpolate(vsf);

    const fvMesh& mesh = ssf.mesh();
    
    // Modify values at the implicit contact boundary
    forAll(mesh.boundary(), patchi)
    {
#ifdef OPENFOAMFOUNDATION        
        // if
        // (
        //     isType<implicitContactFvPatchVectorField>
        //     (
        //         vsf.boundaryField()[patchi]
        //     )
        //     or
        //     isType<implicitGapContactFvPatchVectorField>
        //     (
        //         vsf.boundaryField()[patchi]
        //     )
        // )
        // {
        //     ssf.boundaryFieldRef()[patchi] = vsf.boundaryField()[patchi];
        // }
#endif        
    }

    tmp<GeometricField<GradType, fvPatchField, volMesh>> tgGrad
    (
        gradf(ssf, name)
    );
    GeometricField<GradType, fvPatchField, volMesh>& gGrad = tgGrad.ref();

    correctBoundaryConditions(vsf, gGrad);

    return tgGrad;
}


template<class Type>
void Foam::fv::gaussLinearImplicitContactGrad<Type>::correctBoundaryConditions
(
    const GeometricField<Type, fvPatchField, volMesh>& vsf,
    GeometricField
    <
        typename outerProduct<vector, Type>::type, fvPatchField, volMesh
    >& gGrad
)
{
    typename GeometricField
    <
        typename outerProduct<vector, Type>::type, fvPatchField, volMesh
    >::Boundary& gGradbf = gGrad.boundaryFieldRef();

    forAll(vsf.boundaryField(), patchi)
    {
        if 
        (
            !vsf.boundaryField()[patchi].coupled()
#ifdef OPENFOAMFOUNDATION                    
            // or              
            // isType<implicitContactFvPatchVectorField>
            // (
            //     vsf.boundaryField()[patchi]
            // )
            // or
            // isType<implicitGapContactFvPatchVectorField>
            // (
            //     vsf.boundaryField()[patchi]
            // )
#endif            
        )
        {
            const vectorField n
            (
                vsf.mesh().Sf().boundaryField()[patchi]
              / vsf.mesh().magSf().boundaryField()[patchi]
            );

            gGradbf[patchi] += n *
            (
                vsf.boundaryField()[patchi].snGrad()
              - (n & gGradbf[patchi])
            );
        }
     }
}


// ************************************************************************* //
