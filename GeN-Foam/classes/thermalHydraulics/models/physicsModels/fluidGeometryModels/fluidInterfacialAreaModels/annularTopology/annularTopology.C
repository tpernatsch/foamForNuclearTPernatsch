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

#include "annularTopology.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fluidInterfacialAreaModels
{
    defineTypeNameAndDebug(annularTopology, 0);
    addToRunTimeSelectionTable
    (
        fluidInterfacialAreaModel, 
        annularTopology, 
        fluidInterfacialAreaModels);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidInterfacialAreaModels::annularTopology::annularTopology
(
    const objectRegistry& objReg,
    const dictionary& dict,
    const fluid& dispersed,
    const fluid& continuous
)
:
    fluidInterfacialAreaModel
    (
        objReg,
        dict,
        dispersed,
        continuous
    ),
    DhDispersed_(objReg.lookupObject<volScalarField>("DhDispersed")),
    cutoffAlpha_(dict.lookupOrDefault<scalar>("cutoffAlpha", 0.9))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidInterfacialAreaModels::annularTopology::~annularTopology()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> 
Foam::fluidInterfacialAreaModels::annularTopology::iA() const
{
    tmp<volScalarField> tiA
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("iA", typeName),
                dispersed_.mesh().time().timeName(),
                dispersed_.mesh()
            ),
            dispersed_.mesh(),
            dimensionedScalar("", dimArea/dimVol, 0.0)
        )
    );
    volScalarField& iA = tiA.ref();

    forAll(dispersed_.mesh().cells(), i)
    {
        scalar gamma(dispersed_[i]+continuous_[i]);
        scalar alpha(dispersed_.normalized()[i]);
        scalar& iAi(iA[i]);
        iAi = gamma*4.0*Foam::sqrt(alpha)/max(DhDispersed_[i], 1e-4);
        if (alpha > cutoffAlpha_)
        {
            iAi *= (1.0-alpha)/(1.0-cutoffAlpha_);
        }
    }
    
    return tiA;
}


// ************************************************************************* //
