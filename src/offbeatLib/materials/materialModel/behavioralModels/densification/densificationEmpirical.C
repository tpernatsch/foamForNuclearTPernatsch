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

#include "densificationEmpirical.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(densificationEmpirical, 0);
    addToRunTimeSelectionTable(densificationModel, densificationEmpirical, dictionary);

    const char* densificationEmpirical::group_ = 
        ("behavioralModels::densification::empirical");

    defineParameter(densificationEmpirical, F_epsilonDensification_, 
        "F_epsilonDensification", (dimless), 1.0);
        
    defineParameter(densificationEmpirical, delta_epsilonDensification_, 
        "delta_epsilonDensification", (dimless), 0.0);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::densificationEmpirical::densificationEmpirical
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& baseDict
)
:
    densificationModel(mesh, dict, baseDict),
    burnupName_(dict.lookupOrDefault<word>("burnupName", "Bu" )),
    Bu_(nullptr),
    densityFrac_
    (
        dict.lookupOrDefault<scalar>
        (
            "densityFraction",
            baseDict.lookupOrDefault<scalar>("densityFraction", 0.95)
        )
    ),
    densificationDensityChange_(readScalar(dict.lookup("densificationDensityChange"))),
    densificationTimeConstant_(readScalar(dict.lookup("densificationTimeConstant")))
{    }

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::densificationEmpirical::~densificationEmpirical()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::densificationEmpirical::correct
(
    const scalarField& T, 
    const labelList& addr
)
{
    // Check which F and delta should be used
    const dimensionedScalar& F = 
    isInUserParameters(this->group(), "F_epsilonDensification") ?
    this->F_epsilonDensification() : 
    densificationModel::F_epsilonDensification();

    const dimensionedScalar& delta = 
    isInUserParameters(this->group(), "delta_epsilonDensification") ?
    this->delta_epsilonDensification() : 
    densificationModel::delta_epsilonDensification();

    if(Bu_ == nullptr)
    {        
        Bu_ = &mesh_.lookupObject<volScalarField>(burnupName_);
    }

    const scalarField oldTimeDensification_ = 
    epsilonDensification_.oldTime().internalField();
    
    forAll(addr, addrI)
    {
        const label cellI = addr[addrI];
        
        scalar dLOverL_max(-densificationDensityChange_/(
            densityFrac_*100 + densificationDensityChange_)/3);
        
        // Convert burnup to MWd/tU
        scalar Bui = Bu_->internalField()[cellI]/0.881;
        
        scalar nominalValue =
        dLOverL_max*(1 - exp(-Bui/densificationTimeConstant_)); 
        
        // Perturb epsilon densification for sensitivity analysis
        nominalValue = nominalValue * F.value() + delta.value();     
        
        //- Check if new densification value is lower than previous one           
        //- Keep the smallest value
        scalarField& densificationI = epsilonDensification_.ref();

        densificationI[cellI] = 
        min(nominalValue, oldTimeDensification_[cellI] );
    }
}


// ************************************************************************* //
