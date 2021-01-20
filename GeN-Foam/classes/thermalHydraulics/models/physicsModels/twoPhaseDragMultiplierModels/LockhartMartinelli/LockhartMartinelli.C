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

#include "LockhartMartinelli.H"
#include "fluid.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace twoPhaseDragMultiplierModels
{
    defineTypeNameAndDebug(LockhartMartinelli, 0);
    addToRunTimeSelectionTable
    (
        twoPhaseDragMultiplierModel, 
        LockhartMartinelli, 
        twoPhaseDragMultiplierModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseDragMultiplierModels::LockhartMartinelli::LockhartMartinelli
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
    ),
    C_
    (
        dict.lookupOrDefault<scalar>("C", 20)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::twoPhaseDragMultiplierModels::LockhartMartinelli::~LockhartMartinelli()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::twoPhaseDragMultiplierModels::LockhartMartinelli::correct()
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
        phi2[i] = 1.0 + C_/X[i] + 1.0/sqr(X[i]);
    }

    //- This under-relaxes phi2_ and limits its extrema
    this->setPhi2
    (
        phi2
    );
}


// ************************************************************************* //
