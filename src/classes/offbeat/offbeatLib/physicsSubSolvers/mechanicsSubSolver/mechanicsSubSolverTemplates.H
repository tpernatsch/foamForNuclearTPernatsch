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

#ifndef mechanicsSubSolverTemplates_H
#define mechanicsSubSolverTemplates_H

#include "mechanicsSubSolver.H"
#include "fvc.H"
#include "reverseLinear.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::volTensorField>
Foam::mechanicsSubSolver::correctedSigmaExp
(
    const GeometricField<Type, fvPatchField, volMesh>& sigma,
    const volScalarField& twoMuLambda,
    const volTensorField& gradD
) const
{
    tmp<volTensorField> tsigmaExp
    (
        new volTensorField
        (
            IOobject
            (
                "tsigmaExp",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedTensor
            (
                "tsigmaExp", 
                sigma.dimensions(),
                Foam::tensor::zero
            ),
            // Multi-material derivation imposes a zeroGrad. explicit stress
            multiMaterialCorrection_? 
                zeroGradientFvPatchField<tensor>::typeName
                :
                calculatedFvPatchField<tensor>::typeName
        )
    );

    volTensorField& sigmaExp = tsigmaExp.ref();

    sigmaExp = sigma - twoMuLambda*gradD;

    if(multiMaterialCorrection_)
    {
        // Ensure zeroGradient condition at boundary
        sigmaExp.correctBoundaryConditions();
    }

    return tsigmaExp;
}


template<class Type>
Foam::tmp<Foam::volVectorField>
Foam::mechanicsSubSolver::correctedDivSigmaExp
(
    const volTensorField& sigmaExp,
    const GeometricField<Type, fvPatchField, volMesh>& sigma,
    const volScalarField& twoMuLambda,
    const surfaceScalarField& twoMuLambda_f,
    const volVectorField& D,
    const volTensorField& gradD
) const
{    
    tmp<volVectorField> tdivSigmaExp
    (
        new volVectorField
        (
            IOobject
            (
                "tdivSigmaExp",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector
            (
                "tdivSigmaExp", 
                sigma.dimensions()/dimLength,
                Foam::vector::zero
            )
        )
    );

    volVectorField& divSigmaExp = tdivSigmaExp.ref();

    // Multi-material derivation imposes a RhieChow approach
    if (multiMaterialCorrection_)
    {
        // Normal vectors
        surfaceVectorField nf = mesh().Sf()/mesh().magSf();  

        // Delta vectors
        surfaceVectorField deltaf = mesh().delta();

        // Non-orthogonal correction vectors
        surfaceVectorField k = (nf - deltaf/(nf&deltaf));

        // Reverse-linear interpolated surface gradient 
        surfaceTensorField gradDf = 
            reverseLinear<tensor>(mesh()).interpolate(gradD);

        // Reverse-linear interpolated explicit stress component
        // (the sigmaExp BCs have been set to zeroGradient)
        surfaceTensorField sigma_f = 
        twoMuLambda_f
        *reverseLinear<tensor>(mesh()).interpolate(sigmaExp/twoMuLambda);
        
        divSigmaExp = fvc::div( mesh().Sf() & sigma_f) 
        + fvc::div(twoMuLambda_f*mesh().magSf()* k & gradDf);
    }
    else if (RhieChowCorrection_)
    {
        const scalar scaleFactor
        (
            mechanicsDict_.lookupOrDefault<scalar>("RhieChowScaleFactor", 1.0)
        );

        word laplacianName("laplacian(D" + D.name() + "," + D.name() + ")");

        surfaceTensorField gradDf = 
            fvc::interpolate(gradD, "interpolate(grad("+ D.name() + "))");  

        // Divergence of the explicit component with RhieChow correction.
        // The correction can be partial or total depending on the scale factor
        divSigmaExp = fvc::div(sigma, "div(sigma)") 
        - scaleFactor*fvc::div(twoMuLambda_f*mesh().Sf() & gradDf)
        - (1 - scaleFactor)*fvc::laplacian(twoMuLambda_f, D, laplacianName);
    }
    else
    {
        word laplacianName("laplacian(D" + D.name() + "," + D.name() + ")");

        // Divergence of the explicit component without RhieChow correction
        divSigmaExp = fvc::div(sigma, "div(sigma)") 
            - fvc::laplacian(twoMuLambda_f, D, laplacianName);
    }

    return tdivSigmaExp;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

