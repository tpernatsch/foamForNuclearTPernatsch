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

#if __has_include("commDataLayer.H") 
#include "commDataLayer.H"
#define isCommDataLayerIncluded
#endif

#ifdef isCommDataLayerIncluded


#include "fixedTemperatureFMU.H"
#include "structure.H"
#include "addToRunTimeSelectionTable.H"
#include "commDataLayer.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace powerModels
{
    defineTypeNameAndDebug(fixedTemperatureFMU, 0);
    addToRunTimeSelectionTable
    (
        powerModel, 
        fixedTemperatureFMU, 
        powerModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::powerModels::fixedTemperatureFMU::fixedTemperatureFMU
(
    const structure& structureRef,
    const dictionary& dicts
)
:
    powerModel
    (
        structureRef,
        dicts
    ),
    T_
    (
        IOobject
        (
            "T."+typeName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE //AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimTemperature, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    temperatureNameFromFMU_("temperatureCoupled")
{
    this->setInterfacialArea();
    structure_.setRegionField(*this, T_, "T");


    forAll(this->toc(), regioni)
    {
        word region(this->toc()[regioni]);
        const dictionary& dict(this->subDict(region));

        // Preparing data to update temperature
        word temperatureKeyFromFMU("temperatureNameFromFMU");
        if (dict.found(temperatureKeyFromFMU))
        {
            temperatureNameFromFMU_ = dict.get<word>(temperatureKeyFromFMU);
            // Communicating with the FMU
            const Time& runTime = this->db().time();
            commDataLayer& data = commDataLayer::New(runTime); 
            // Store in data layer and set its initial value to the T 
            // in the dictionary      
            data.storeObj(
                dict.get<scalar>("T"),
                temperatureNameFromFMU_,
                commDataLayer::causality::in
                );
            Info << "Using FMUs for the temperature in " << dict.dictName() << endl;
        }
    }

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::powerModels::fixedTemperatureFMU::~fixedTemperatureFMU()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::powerModels::fixedTemperatureFMU::temperatureUpdate() const
{
    forAll(this->toc(), regioni)
    {
        // Communicating with the FMU
        const Time& runTime = this->db().time();
        commDataLayer& data = commDataLayer::New(runTime);
        const scalar temperatureFromFMU =
            data.getObj<scalar>(temperatureNameFromFMU_,commDataLayer::causality::in);
        
        //- Setup cellToRegion_ mapping
        word region(this->toc()[regioni]);
        const labelList& regionCells
        (
            structure_.cellLists()[region]
        );

        forAll(regionCells, i)
        {
            label celli(regionCells[i]);
            T_[celli] =  temperatureFromFMU;
        }          

    }
}

void Foam::powerModels::fixedTemperatureFMU::correctT(volScalarField& T) const
{
    this->temperatureUpdate();
    forAll(cellList_, i)
    {
        label celli(cellList_[i]);
        T[celli] = T_[celli];
    }
}

void Foam::powerModels::fixedTemperatureFMU::powerOff()
{
    //- If you set iA to 0, the energy contribution from this powerModel to the
    //  fluid energy equation will be 0, equivalent to a "power" off scenario
    iA_ *= 0.0;
}

#endif
// ************************************************************************* //
