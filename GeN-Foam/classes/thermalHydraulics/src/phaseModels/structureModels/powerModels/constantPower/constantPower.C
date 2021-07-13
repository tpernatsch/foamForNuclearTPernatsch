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

#include "constantPower.H"
#include "addToRunTimeSelectionTable.H"
#include "structure.H"
#include "fvMatrix.H"
#include "fvmDdt.H"
#include "fvmSup.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace powerModels
{
    defineTypeNameAndDebug(constantPower, 0);
    addToRunTimeSelectionTable
    (
        powerModel, 
        constantPower, 
        powerModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::powerModels::constantPower::constantPower
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
    powerDensity_
    (
        IOobject
        (
            IOobject::groupName("powerDensity", typeName),
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimPower/dimVol, 0.0),
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
    ),
    alphaRhoCp_
    (
        IOobject
        (
            IOobject::groupName("alphaRhoCp", typeName),
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("", dimEnergy/dimVol/dimTemperature, 0.0),
        zeroGradientFvPatchScalarField::typeName
    )
{
    this->setInterfacialArea();
    
    structure_.setRegionField(this, powerDensity_, "powerDensity");
    structure_.setRegionField(this, T_, "T");

    forAll(this->toc(), regioni)
    {
        word region(this->toc()[regioni]);
        const dictionary& dict(this->subDict(region));
        
        //- Setup cellToRegion_ mapping
        const labelList& regionCells
        (
            structure_.cellLists()[region]
        );
        scalar rhoCp(0);
        if (dict.found("rho") and dict.found("Cp"))
        {
            rhoCp =
                dict.get<scalar>("rho")*
                dict.get<scalar>("Cp");
        }
        else if (dict.found("rhoCp"))
        {
            rhoCp = dict.get<scalar>("rhoCp");
        }
        else
        {
            FatalErrorInFunction
                << "constantPower region: " << region << " -> "
                << "specify either rhoCp or both rho and Cp"
                << exit(FatalError);
        }
        forAll(regionCells, i)
        {
            label celli(regionCells[i]);
            alphaRhoCp_[celli] = rhoCp;
        }
    }

    //- The alphaRhoCp is read as a rhoCp, alpha is multiplied at this step
    alphaRhoCp_ = 
        max
        (
            alpha_*alphaRhoCp_, 
            dimensionedScalar
            ("", dimEnergy/dimVol/dimTemperature, 1e-69)
        );
    alphaRhoCp_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::powerModels::constantPower::~constantPower()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::powerModels::constantPower::correctT(volScalarField& T) const
{
    forAll(cellList_, i)
    {
        label celli(cellList_[i]);
        T[celli] = T_[celli];
    }
    //- Correct boundary conditions called on T in structure.C after it has
    //  been set, cell-by-cell, by all powerModels
}


void Foam::powerModels::constantPower::powerOff()
{
    powerDensity_ *= 0.0;
}


void Foam::powerModels::constantPower::correct
(
    const volScalarField& HTSum,  // == SUM_j [htc_j*T_j*frac_j]
    const volScalarField& HSum    // == SUM_j [htc_j*frac_j]
)
{
    scalar dt(mesh_.time().deltaT().value());
    volScalarField& T0(T_.oldTime()); 
    
    forAll(cellList_, i)
    {
        label celli(cellList_[i]);
        scalar alphaRhoCpByDt(alphaRhoCp_[celli]/dt);
        const scalar& iA(iA_[celli]);
        T_[celli] = 
            (
                iA*HTSum[celli] 
            +   alpha_[celli]*powerDensity_[celli]
            +   alphaRhoCpByDt*T0[celli]
            )/
            (alphaRhoCpByDt + iA*HSum[celli]);
    }
    T_.correctBoundaryConditions();

    /*
    What I did above is equivalent to doing this, but possibly faster, given
    that I might not have to solve the equation for the entire mesh, but only
    over the cells over which the constantPower models exists. The time
    derivative was hardcoded to be forward Euler

    fvScalarMatrix EEqn
    (
            fvm::ddt(alphaRhoCp_, T_)
        ==
            iA_*HTSum
        -   fvm::Sp(iA_*HSum, T_)
        +   alpha_*powerDensity_
    );
    EEqn.solve();
    T_.correctBoundaryConditions();
    */
}


// ************************************************************************* //
