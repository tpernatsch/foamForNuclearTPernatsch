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
#include "zeroGradientFvPatchFields.H"
#include "fluid.H"
#include "structureModel.H"
#include "fvCFD.H"

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
    const volScalarFieldPtrTable& htcs,
    volScalarField& dmdt,
    volScalarField& iT,
    volScalarField& iA
)
:
    IOdictionary
    (
        IOobject
        (
            IOobject::groupName("phaseChangeModel", typeName),
            fluid1.mesh().time().constant(),
            fluid1.mesh()
        ),
        dict
    ),
    mesh_(fluid1.mesh()),
    pimple_(pimple),
    fluid1_(fluid1),
    fluid2_(fluid2),
    structure_(mesh_.lookupObject<structureModel>("alpha.structure")),
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
            1e-6
        )
    ),
    residualIACells_(0)
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

    //- Knowledge of which phase is liquid and which is vapour are not
    //  required by all phaseChange models (e.g. the heatDrivenPhaseChange one
    //  does not require it and distinguishes the phases precisely on the basis
    //  of the latent heat (pardon the inconsistency with the error message,
    //  yet it's for the user's sake). Having this check regardless of the
    //  models adds however a further layer of double-checking thinks that can
    //  be beneficial
    if 
    (
        fluid1_.isLiquid() and fluid2_.isGas() and
        (
            fvc::domainIntegrate(fluid1_.thermo().hc()).value() >=
            fvc::domainIntegrate(fluid2_.thermo().hc()).value() 
        )
    )
    {
        FatalErrorInFunction
            << "Negative latent heat! Ensure that the liquid phase (i.e. " 
            << fluid1_.name() << ") enthalpy of formation (Hc) in its "
            << "thermoPhysicalProperties file is smaller than the vapour "
            << "phase (i.e. " << fluid2_.name() << ") enthalpy of formation" 
            << exit(FatalError);
    }
    else if 
    (
        fluid1_.isGas() and fluid2_.isLiquid() and
        (
            fvc::domainIntegrate(fluid1_.thermo().hc()).value() <=
            fvc::domainIntegrate(fluid2_.thermo().hc()).value() 
        )
    )
    {
        FatalErrorInFunction
            << "Negative latent heat! Ensure that the liquid phase (i.e. " 
            << fluid2_.name() << ") enthalpy of formation (Hc) in its "
            << "thermoPhysicalProperties file is smaller than the vapour "
            << "phase (i.e. " << fluid1_.name() << ") enthalpy of formation" 
            << exit(FatalError);
    }
    else if 
    (
        (!fluid1_.isLiquid() and !fluid1_.isGas()) or 
        (!fluid2_.isLiquid() and !fluid2_.isGas())
    )
    {
        FatalErrorInFunction
            << "Either phase " << fluid1_.name() << " or " << fluid2_.name()
            << " have an undetermined stateOfMatter (should be specified in "
            << "phaseProperties." << fluid1_.name() << "Properties and/or "
            << "phaseProperties." << fluid2_.name() << "Properties)"
            << exit(FatalError);
    }

    //- Read residualIACells from regions
    wordList residualIARegions
    (
        this->lookupOrDefault<wordList>("residualInterfacialAreaRegions", wordList())
    );
    forAll(residualIARegions, i)
    {
        const labelList& regionCells
        (
            structure_.cellLists()[residualIARegions[i]]
        );

        forAll(regionCells, j)
        {
            residualIACells_.append(regionCells[j]);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::phaseChangeModel::relaxationFactor() const
{
    //- Read relaxation factor if present
    scalar f(1.0);
    if (mesh_.relaxField("massTransfer"))
    {
        f = mesh_.fieldRelaxationFactor("massTransfer");
        //- Basically stabilizedMassTransfer amounts to a constant
        //  underrelaxation, even between consecutive time steps and on the 
        //  last PIMPLE iteration (which is not how underrelaxation should be 
        //  applied). However, if the time-steps are small enough (and they do 
        //  get so if the boiling is violent enough, i.e. ~ 10-100  
        //  microseconds), it greatly enchances stability (if the chosen 
        //  relaxation factor is small enough, e.g. ~ 0.1-0.2) and does not
        //  appear to affect results at all.
        bool stabilizedMassTransfer
        (
            pimple_.dict().lookupOrDefault<bool>
            (
                "stabilizedMassTransfer", false
            )
        );

        if 
        (
            !stabilizedMassTransfer
        and (pimple_.firstIter() or pimple_.finalIter())
        )   f = 1.0;
    }

    return f;
}


void Foam::phaseChangeModel::limitInterfacialArea()
{
    if (residualIA_.value() != 0.0)
    {
        if (residualIACells_.size() != 0)
        {
            scalar residualIAValue(residualIA_.value());
            forAll(residualIACells_, i)
            {
                scalar& iAi(iA_[residualIACells_[i]]);
                iAi = max(iAi, residualIAValue);
            }
        }
        else
        {
            iA_ = max(iA_, residualIA_);
        }
    }
}


void Foam::phaseChangeModel::limitMassTransfer()
{
    //- Prevent removing mass from phases that are not present in a cell
    dmdt_ = posPart(dmdt_)*pos(fluid1_) + negPart(dmdt_)*pos(fluid2_);
}


void Foam::phaseChangeModel::correctInterfacialT()
{
    scalar f(1.0);
    if (mesh_.relaxField("interfacialTemperature"))
    {
        f = mesh_.fieldRelaxationFactor("interfacialTemperature");
    }
    iT_ = (1-f)*iT_ + f*saturation_->Tsat(p_);
}


// ************************************************************************* //
