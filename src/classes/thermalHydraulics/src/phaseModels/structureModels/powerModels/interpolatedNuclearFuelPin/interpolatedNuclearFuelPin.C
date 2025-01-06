/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2412                                                  |
|    Copyright 2011-2016 OpenFOAM Foundation, 2017-2024 OpenCFD Ltd.          |
-------------------------------------------------------------------------------
License
    This file is part of GeN-Foam.

    GeN-Foam is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    GeN-Foam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    This offering is not approved or endorsed by the OpenFOAM Foundation nor
    OpenCFD Limited, producer and distributor of the OpenFOAM(R)software via
    www.openfoam.com, and owner of the OPENFOAM(R) and OpenCFD(R) trademarks.

    This particular snippet of code is developed according to the developer's
    knowledge and experience in OpenFOAM. The users should be aware that
    there is a chance of bugs in the code, though we've thoroughly test it.
    The source code may not be in the OpenFOAM coding style, and it might not
    be making use of inheritance of classes to full extent.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "interpolatedNuclearFuelPin.H"
#include "structure.H"
#include "addToRunTimeSelectionTable.H"
#include "SquareMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace powerModels
{
    defineTypeNameAndDebug(interpolatedNuclearFuelPin, 0);
    addToRunTimeSelectionTable
    (
        powerModel,
        interpolatedNuclearFuelPin,
        powerModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::powerModels::interpolatedNuclearFuelPin::interpolatedNuclearFuelPin
(
    structure& structureRef,
    const dictionary& dicts
)
:
    powerModel
    (
        structureRef,
        dicts
    ),
    rfo_(0),//this->get<scalar>("fuelOuterRadius")),
    cellToRegion_(mesh_.cells().size(), 0),
    regionIndexToRegionName_(0),
    interpolatedHeatFlux_
    (
        IOobject
        (
            "interpolatedHeatFlux",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimPower/dimArea, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    interpolatedPowerDensity_
    (
        IOobject
        (
            "interpolatedPowerDensity",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        interpolatedHeatFlux_*iA_,
        zeroGradientFvPatchScalarField::typeName
    ),
    T_
    (
        IOobject
        (
            IOobject::groupName("T", typeName),
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimTemperature, 0.0),
        zeroGradientFvPatchScalarField::typeName
    )
{

    structure_.setRegionField(*this, interpolatedPowerDensity_, "powerDensity");
    structure_.setRegionField(*this, T_, "T");


    forAll(this->toc(), regioni)
    {
        word region(this->toc()[regioni]);
        const dictionary& dict(this->subDict(region));

        //- Setup cellToRegion_ mapping
        const labelList& regionCells
        (
            structure_.cellLists()[region]
        );
        forAll(regionCells, i)
        {
            label celli(regionCells[i]);
            cellToRegion_[celli] = regioni;
        }

        //- Add to regionIndexToRegionName_ mapping
        regionIndexToRegionName_.append(region);

        //- Read region dict entries
        scalar rfo(dict.get<scalar>("fuelOuterRadius"));
        rfo_.append(rfo);
    }
    this->setInterfacialArea();

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::powerModels::interpolatedNuclearFuelPin::~interpolatedNuclearFuelPin()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::powerModels::interpolatedNuclearFuelPin::setInterfacialArea()
{
    forAll(this->cellList_, i)
    {
        const label& celli(this->cellList_[i]);
        iA_[celli] = 2.0*alpha_[celli]/(rfo_[cellToRegion_[celli]]);
    }
    iA_.correctBoundaryConditions();
}


void Foam::powerModels::interpolatedNuclearFuelPin::correct
(
    const volScalarField& HTSum,  // == SUM_j [htc_j*T_j*frac_j]
    const volScalarField& HSum    // == SUM_j [htc_j*frac_j]
)
{

    interpolatedPowerDensity_=interpolatedHeatFlux_*iA_/alpha_;
    interpolatedPowerDensity_.correctBoundaryConditions();
    forAll(this->cellList_, i)
    {
        const label celli(this->cellList_[i]);
        T_[celli]=(interpolatedHeatFlux_[celli]+HTSum[celli])/HSum[celli];

    }

    T_.correctBoundaryConditions();

}


void Foam::powerModels::interpolatedNuclearFuelPin::correctT(volScalarField& T) const
{
    //- Set T to pin surface temperature, i.e. Tco_
    forAll(cellList_, i)
    {
        label celli(cellList_[i]);
        T[celli] = T_[celli];
    }
}


// ************************************************************************* //
