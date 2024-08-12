/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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

#include "eig3Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

void eig3Field
(
    const volTensorField& A, volTensorField& V, volVectorField& d
)
{
    const tensorField& AI = A.internalField();
    tensorField& VI = V.primitiveFieldRef();
    vectorField& dI = d.primitiveFieldRef();

    forAll(AI, cellI)
    {
        eig3().eigen_decomposition(AI[cellI], VI[cellI], dI[cellI]);
    }

    forAll(A.boundaryField(), patchI)
    {
        if (A.boundaryField()[patchI].type() != "empty")
        {
            const tensorField& AB = A.boundaryField()[patchI];
            tensorField& VB = V.boundaryFieldRef()[patchI];
            vectorField& dB = d.boundaryFieldRef()[patchI];

            forAll(AB, faceI)
            {
                eig3().eigen_decomposition(AB[faceI], VB[faceI], dB[faceI]);
            }
        }
    }
}


void eig3Field
(
    const volSymmTensorField& A, volTensorField& V, volVectorField& d
)
{
    const symmTensorField& AI = A.internalField();
    tensorField& VI = V.primitiveFieldRef();
    vectorField& dI = d.primitiveFieldRef();

    forAll(AI, cellI)
    {
        eig3().eigen_decomposition(AI[cellI], VI[cellI], dI[cellI]);
    }

    forAll(A.boundaryField(), patchI)
    {
        if (A.boundaryField()[patchI].type() != "empty")
        {
            const symmTensorField& AB = A.boundaryField()[patchI];
            tensorField& VB = V.boundaryFieldRef()[patchI];
            vectorField& dB = d.boundaryFieldRef()[patchI];

            forAll(AB, faceI)
            {
                eig3().eigen_decomposition(AB[faceI], VB[faceI], dB[faceI]);
            }
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
