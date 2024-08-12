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

#include "fuelMaterial.H"
#include "addToRunTimeSelectionTable.H"
#include "globalFieldLists.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fuelMaterial, 0);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fuelMaterial::fuelMaterial
(
    const fvMesh& mesh, 
    const dictionary& materialModelDict,
    const labelList& addr
)
:
    materialModel(mesh, materialModelDict, addr),
    densification_(),
    swelling_(),
    relocation_(),
    failure_(),
    poreVelocity_(),
    enrichment_(readScalar(materialModelDict_.lookup("enrichment"))),
    rGrain_(readScalar(materialModelDict_.lookup("rGrain"))),
    oxygenMetalRatio_(materialModelDict_.lookupOrDefault("oxygenMetalRatio", 2.0)),
    densityFrac_(readScalar(materialModelDict_.lookup("densityFraction"))),
    porosity_
    (
        createOrLookup<scalar>
        (
            mesh, 
            "porosity", 
            dimless, 
            0.0,
            "zeroGradient"
        )
    ),
    isotopesU_
    (
        materialModelDict_.lookupOrDefault<scalarList>
        (
            "isotopesU", scalarList(4,0.0)
        )
    ),
    isotopesPu_
    (
        materialModelDict_.lookupOrDefault<scalarList>
        (
            "isotopesPu", scalarList(5,0.0)
        )
    ),
    isotopesAm_
    (
        materialModelDict_.lookupOrDefault<scalarList>
        (
            "isotopesAm", scalarList(2,0.0)
        )
    ),
    ratioUMetal_(materialModelDict_.lookupOrDefault("ratioUMetal", 0.0)),
    ratioPuMetal_(materialModelDict_.lookupOrDefault("ratioPuMetal", 0.0)),
    ratioAmMetal_(materialModelDict_.lookupOrDefault("ratioAmMetal", 0.0)),
    dishFraction_(readScalar(materialModelDict_.lookup("dishFraction", 0.0)))
{
    // Initialize porosity
    forAll(addr_, i)
    {
        const label cellI = addr_[i];   

        // Initialize only if not already set
        if(porosity_[cellI] < VSMALL)
        {
            porosity_[cellI] = 1 - densityFrac_;
        }
    }
 
    porosity_.correctBoundaryConditions();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fuelMaterial::~fuelMaterial()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// ************************************************************************* //
