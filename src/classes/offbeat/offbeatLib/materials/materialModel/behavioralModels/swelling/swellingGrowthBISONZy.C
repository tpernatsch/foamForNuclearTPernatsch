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

#include "swellingGrowthBISONZy.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(swellingGrowthBISONZy, 0);
    addToRunTimeSelectionTable(swellingModel, swellingGrowthBISONZy, dictionary);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::swellingGrowthBISONZy::swellingGrowthBISONZy
(
    const fvMesh& mesh,
    const dictionary& dict 
)
:
    swellingModel(mesh, dict),
    fastFluenceName_(dict.lookupOrDefault<word>("fastFluenceName", "fastFluence" )),
    fastFluence_(nullptr),
    A(3e-20),
    n(0.794),
    cladType_("ESCORE"),
    perturb(1.0)
{    
    if(dict.found("swelling"))
    {
        const dictionary& swellingDict = dict.subDict("swelling");

        A = swellingDict.lookupOrDefault<scalar>("A", 3e-20);
        n = swellingDict.lookupOrDefault<scalar>("n", 0.794);

        perturb = swellingDict.lookupOrDefault<scalar>("perturb", 1.0);
    }  

    dictionary swellingDict
    (
        dict.subOrEmptyDict("swellingModelOptions")
    );

    if(swellingDict.found("cladType"))
    {
        cladType_ = word(swellingDict.lookup("cladType"));

        if (cladType_ == "SRA")
        {
            A = 2.18e-21;
            n = 0.845;
        }
        else if (cladType_ == "RXA")
        {
            A = 1.09e-21;
            n = 0.845;
        }
        else if (cladType_ == "PRA")
        {
            A = 1.09e-21;
            n = 0.845;
        }
        else if (cladType_ == "ZIRLO")
        {
            A = 9.7893e-25;
            n = 0.98239;
        }
        else if (cladType_ == "ESCORE")
        {
            A = 3.0e-20;
            n = 0.794;
        }
        else if (cladType_ == "M5")
        {
            A = 7.013e-21;
            n = 0.81787;
        }
        else
        {
            FatalErrorInFunction
                << nl
                << "    Wrong \"cladType\" " << cladType_ << " for LimbackCreepModel."
                << "Available choices are:" << nl
                << "    - SRA i.e. stress relief annealed (Zr2 or Zr4), default choice" << nl
                << "    - RXA i.e. recrystallization annealed (Zr2 or M5)" << nl
                << "    - PRA i.e. partially recrystallization annealed (Zr2)" << nl
                << "    - ZIRLO i.e. stress relief annealed ZIRLO" << nl
                << "    - ESCORE i.e. model from Rashid" << nl
                << "    - M5 i.e. model from Gilbon" << nl
                << abort(FatalError);
        }
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::swellingGrowthBISONZy::~swellingGrowthBISONZy()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::swellingGrowthBISONZy::correct
(
    const scalarField& T, 
    const labelList& addr
)
{
    if(fastFluence_ == nullptr)
    {        
        fastFluence_ = &mesh_.lookupObject<volScalarField>(fastFluenceName_);
    }

    symmTensorField& swellingI = epsilonSwelling_.ref();
    const symmTensorField& swellingIOld = epsilonSwelling_.oldTime().internalField();
    
    const scalarField& fastFluenceI = 
    fastFluence_->internalField();
    const scalarField& fastFluenceIOld = 
    fastFluence_->oldTime().internalField();

    forAll(addr, addrI)
    {
        const label cellI = addr[addrI];

        const scalar fastFluence = fastFluenceI[cellI];
        const scalar fastFluenceOld = fastFluenceIOld[cellI];
        
        const scalar nominalValue = A*pow(fastFluence,n);  
        const scalar nominalValueOld = A*pow(fastFluenceOld,n);  
        const scalar swellIncrement = nominalValue - nominalValueOld;

        swellingI[cellI] = swellingIOld[cellI];

        swellingI[cellI].xx() += -(1 - pow((1 + swellIncrement*perturb), -0.5));  
        swellingI[cellI].yy() += -(1 - pow((1 + swellIncrement*perturb), -0.5));
        swellingI[cellI].zz() += (swellIncrement*perturb);  
    }

    // Check which F and delta should be used
    const dimensionedScalar& F = 
    isInUserParameters(this->group(), "F_epsilonSwelling") ?
    this->F_epsilonSwelling() : 
    swellingModel::F_epsilonSwelling();

    const dimensionedScalar& delta = 
    isInUserParameters(this->group(), "delta_epsilonSwelling") ?
    this->delta_epsilonSwelling() : 
    swellingModel::delta_epsilonSwelling();

    // Perturb epsilon swelling for sensitivity analysis
    applyParameters(epsilonSwelling_.ref(), addr, F, delta);

}


// ************************************************************************* //
