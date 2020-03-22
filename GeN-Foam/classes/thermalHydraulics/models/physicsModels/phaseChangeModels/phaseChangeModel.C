/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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

#include "phaseChangeModel.H"
#include "fluid.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseChangeModel, 0);
    defineRunTimeSelectionTable(phaseChangeModel, phaseChangeModels);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeModel::phaseChangeModel
(
    const dictionary& dict,
    const pimpleControl& pimple, 
    const fluid& fluid1,
    const fluid& fluid2,
    const volScalarField& p,
    const volScalarFieldTable& htcs,
    volScalarField& dmdt,
    volScalarField& iT,
    volScalarField& iA
)
:
    IOdictionary
    (
        IOobject
        (
            "phaseChangeModel",
            fluid1.mesh().time().constant(),
            fluid1.mesh()
        ),
        dict
    ),
    mesh_(fluid1.mesh()),
    pimple_(pimple),
    fluid1_(fluid1),
    fluid2_(fluid2),
    p_(p),
    htcs_(htcs),
    dmdt_(dmdt),
    iT_(iT),
    iA_(iA),
    saturation_
    (
        saturationModel::New
        (
            this->subDict("saturationModel"),
            fluid1_,
            fluid2_
        )
    ),
    relaxOnFinalIter_
    (
        this->lookupOrDefault<bool>("relaxOnFinalIter", false)
    ),
    residualIA_
    (
        dimensionedScalar::lookupOrDefault
        (
            "residualInterfacialArea",
            *this,
            dimArea/dimVolume,
            1e-9
        )
    )
{
    //- Set initial interfacial temperature (to avoid problems with the
    //  under-relaxation at the first time-step in EEqns.H
    IOobject iTHeader
    (
        iT_.name(),
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ
    );
    if (!iTHeader.typeHeaderOk<volScalarField>(true))
    {
        iT_ = saturation_->Tsat(p_);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::phaseChangeModel::relaxationFactor() const
{
    scalar f(mesh_.fieldRelaxationFactor("dmdt"));
    if (pimple_.finalIter() and !relaxOnFinalIter_)
    {
        f = 1.0;
    }

    return f;
}


void Foam::phaseChangeModel::limitInterfacialArea()
{
    if (residualIA_.value() != 0.0)
    {
        iA_ = max(iA_, residualIA_);
    }
}


void Foam::phaseChangeModel::limitMassTransfer()
{
    //- Prevent removing mass from phases that are not present in a cell
    dmdt_ = posPart(dmdt_)*pos(fluid1_) + negPart(dmdt_)*pos(fluid2_);
}


void Foam::phaseChangeModel::correctInterfacialT()
{
    scalar f(mesh_.fieldRelaxationFactor("iT12"));
    iT_ = (1-f)*iT_ + f*saturation_->Tsat(p_);
}


// ************************************************************************* //
