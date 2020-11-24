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

#include "KottowskiSavatteri.H"
#include "fluid.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace twoPhaseDragMultiplierModels
{
    defineTypeNameAndDebug(KottowskiSavatteri, 0);
    addToRunTimeSelectionTable
    (
        twoPhaseDragMultiplierModel, 
        KottowskiSavatteri, 
        twoPhaseDragMultiplierModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseDragMultiplierModels::KottowskiSavatteri::KottowskiSavatteri
(
    const objectRegistry& objReg,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    twoPhaseDragMultiplierModel
    (
        objReg,
        dict,
        mesh
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::twoPhaseDragMultiplierModels::KottowskiSavatteri::~KottowskiSavatteri()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::twoPhaseDragMultiplierModels::KottowskiSavatteri::correct()
{
    volScalarField phi2
    (
        IOobject
        (
            "",
            mesh_.time().timeName(),
            this->db()
        ),
        mesh_,
        dimensionedScalar("", dimless, 1)
    );

    const volScalarField& X(mFluid_.XLM());
    forAll(X, i)
    {
        //- The correlation is valid only for 7e-2 < X < 30
        scalar log10X(log10(min(max(X[i], 0.07), 30)));
        phi2[i] = 
            pow
            (
                10,
                2.0*
                (
                    0.1046*sqr(log10X)
                -   0.5098*log10X
                +   0.6252
                )
            );
    }

    //- This under-relaxes phi2_ and limits its extrema
    this->setPhi2
    (
        phi2
    );
}


// ************************************************************************* //
