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

#include "constantSwelling.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(constantSwelling, 0);
    addToRunTimeSelectionTable(swellingModel, constantSwelling, dictionary);

    const char* constantSwelling::group_ = 
        ("behavioralModels::swelling::constant");

    defineParameter(constantSwelling, F_epsilonSwelling_, 
        "F_epsilonSwelling", (dimless), 1.0);
        
    defineParameter(constantSwelling, delta_epsilonSwelling_, 
        "delta_epsilonSwelling", (dimless), 0.0);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constantSwelling::constantSwelling
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& baseDict
)
:
    swellingModel(mesh, dict, baseDict),
    s_(readScalar(dict.lookup("swellingRate"))),
    fastFluenceName_(dict.lookupOrDefault<word>("fastFluenceName", "fastFluence" )),
    fastFluence_(nullptr)
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::constantSwelling::~constantSwelling()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::constantSwelling::correct
(
    const scalarField& T,
    const labelList& addr
)
{
    // Check which F and delta should be used
    const dimensionedScalar& F =
    isInUserParameters(this->group(), "F_epsilonSwelling") ?
    this->F_epsilonSwelling() :
    swellingModel::F_epsilonSwelling();

    const dimensionedScalar& delta =
    isInUserParameters(this->group(), "delta_epsilonSwelling") ?
    this->delta_epsilonSwelling() :
    swellingModel::delta_epsilonSwelling();

    const fvMesh& mesh_ = epsilonSwelling_.mesh();

    if(fastFluence_ == nullptr)
    {
        fastFluence_ = &mesh_.lookupObject<volScalarField>(fastFluenceName_);
    }

    symmTensorField& swellingI = epsilonSwelling_.ref();

    const volScalarField& phii = *fastFluence_;

    forAll(addr, addrI)
    {
        const label cellI = addr[addrI];
        //- Unit at  10e25 n/m^2
        const scalar phi = phii[cellI]*1e4/1e25;

        scalar nominalValue(s_*phi);

        //Perturb for sensitivity analysis
        nominalValue = nominalValue * F.value() + delta.value();

        swellingI[cellI] = nominalValue * I;
    }


}


// ************************************************************************* //
