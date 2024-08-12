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

#include "burnupFromPower.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(burnupFromPower, 0);
    addToRunTimeSelectionTable
    (
        burnup, 
        burnupFromPower, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::burnupFromPower::burnupFromPower
(
    const fvMesh& mesh,
    const materials& mat,
    const dictionary& burnupDict
)
:
    constantBurnup(mesh, mat, burnupDict),
    heatSourceName_("Q")
{  
    if(      
#ifdef OPENFOAMFOUNDATION
        !burnupDict.isNull()
#elif OPENFOAMESI
        !burnupDict.isNullDict()
#endif  
    )      
    {
        heatSourceName_ = 
        burnupDict.lookupOrDefault<word>("heatSourceName", "Q");
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::burnupFromPower::~burnupFromPower()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::burnupFromPower::correct()
{   
    // TODO: Check whether the Bu update should be done at each iteration or
    // only at the beginning of the time step
    const volScalarField& Q_ = 
    mesh_.lookupObject<volScalarField>(heatSourceName_);

    const scalarField& rhoi = mat_.rho().internalField();
    const scalarField& Qi = Q_.internalField();
    const scalarField& BuOldi = Bu_.oldTime().internalField();
    scalarField& Bui = Bu_.ref();        

    // Dt transformed to be in days
    scalar deltaT = mesh_.time().deltaT().value()/3600.0/24.0;

    forAll(Qi, i)
    {
        // Bu in MWd/MTUO2   
        Bui[i] = BuOldi[i] + (Qi[i]*deltaT/rhoi[i]/1000);
    }
    
    Bu_.correctBoundaryConditions();
}
// ************************************************************************* //
