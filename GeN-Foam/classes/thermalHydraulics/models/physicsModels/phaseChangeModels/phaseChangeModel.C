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
    latentHeat_
    (
        latentHeatModel::New
        (
            this->subDict("latentHeatModel"),
            fluid1_,
            fluid2_,
            p,
            dmdt,
            iT
        )
    ),
    saturation_
    (
        saturationModel::New
        (
            this->subDict("saturationModel"),
            fluid1_,
            fluid2_
        )
    ),
    residualIA_
    (
        dimensionedScalar::lookupOrDefault
        (
            "residualInterfacialArea",
            *this,
            dimArea/dimVolume,
            1e-3
        )
    ),
    residualIACells_(0),
    dmLostToLimiter_(0),
    posDmdtResidualIA_
    (
        dimensionedScalar::lookupOrDefault
        (
            "posDmdtResidualInterfacialArea",
            *this,
            dimArea/dimVolume,
            0.0
        )
    ),
    negDmdtResidualIA_
    (
        dimensionedScalar::lookupOrDefault
        (
            "negDmdtResidualInterfacialArea",
            *this,
            dimArea/dimVolume,
            0.0
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

    //- 
    if 
    (
        !(fluid1_.isLiquid() and fluid2_.isGas())
    and !(fluid2_.isLiquid() and fluid1_.isGas())
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
            if (posDmdtResidualIA_.value() != 0.0)
            {
                forAll(dmdt_, i)
                {
                    scalar& iA(iA_[i]);
                    if (dmdt_[i] > 0.0)
                        iA = max(iA, posDmdtResidualIA_.value());
                }
            }
            if (negDmdtResidualIA_.value() != 0.0)
            {
                forAll(dmdt_, i)
                {
                    scalar& iA(iA_[i]);
                    if (dmdt_[i] < 0.0)
                        iA = max(iA, negDmdtResidualIA_.value());
                }
            }
        }
    }
}


void Foam::phaseChangeModel::limitMassTransfer()
{
    //- Prevent removing mass from phases that are not present in a cell
    dmdt_ = posPart(dmdt_)*pos(fluid1_) + negPart(dmdt_)*pos(fluid2_);

    /*
    Given that the continuity equations are:

    ddt(alpha1*rho1) + div(alphaRhoPhi1) = -dmdt
    ddt(alpha2*rho2) + div(alphaRhoPhi2) = dmdt

    and that they are solved explicitly (albeit with the whole MULES limiter
    step before the explicit solution), there should be a maximum allowable
    limit on the value of the mass transfer term. This can be computed from 
    the continuity equations assuming the maximum allowable ddt within a 
    single time-step (i.e. that ddt that would make my phase fraction change
    from its current value to 0 in one time step). Due to the way the 
    equations are implemented, a maximum on dmdt is derived from the 
    continuity equation of phase 1 while a minimum is derived from that of
    phase 2. Continuity errors are (for now) not considered.
    */

    if (pimple_.dict().found("massTransferSafetyFactor"))
    {
        scalar SF(pimple_.dict().get<scalar>("massTransferSafetyFactor"));

        scalar dt(mesh_.time().deltaT().value());
        const scalarField& V(mesh_.V());
        scalar totV(0);
        forAll(V, i)
        {
            totV += V[i];
        }

        scalarField div1 = fvc::div(fluid1_.alphaRhoPhi())().primitiveField();
        scalarField div2 = fvc::div(fluid2_.alphaRhoPhi())().primitiveField();

        forAll(dmdt_, i)
        {
            scalar& dmdt(dmdt_[i]);
            if (dmdt > 0.0)
            {
                const scalar& Vi(V[i]);
                scalar maxDmdt
                (
                    max(SF*fluid1_[i]*fluid1_.rho()[i]/dt-div1[i], 0.0)
                );
                scalar dmdt0(dmdt);
                dmdt = min(dmdt, maxDmdt);
                if (pimple_.finalIter())
                    dmLostToLimiter_ += (dmdt0-dmdt)*Vi*dt;
            }
            else if (dmdt < 0.0)
            {
                const scalar& Vi(V[i]);
                scalar minDmdt
                (
                    min(-SF*fluid2_[i]*fluid2_.rho()[i]/dt+div2[i], 0.0)
                );
                scalar dmdt0(dmdt);
                dmdt = max(dmdt, minDmdt);
                if (pimple_.finalIter())
                    dmLostToLimiter_ += (dmdt0-dmdt)*Vi*dt;
            }
        }

        if (pimple_.finalIter())
        {
            Info<< "Cumulative dm lost to limiter = " 
                << (dmLostToLimiter_) 
                << " kg/m3" << endl;
        }
    }

    dmdt_.correctBoundaryConditions();
}


void Foam::phaseChangeModel::correctInterfacialT()
{
    if (mesh_.relaxField("T.interface"))
        iT_.storePrevIter();

    if (pimple_.dict().found("maxTInterfaceRelChange"))
    {
        scalar maxRelChange
        (
            pimple_.dict().get<scalar>
            (
                "maxTInterfaceRelChange"
            )
        );
        scalarField relChange
        (
            min
            (
                max
                (
                    saturation_->Tsat(p_)().primitiveField()/
                    iT_.primitiveField(),
                    1.0 - maxRelChange
                ),
                1.0 + maxRelChange
            )
        );
        iT_.primitiveFieldRef() *= relChange;
        iT_.correctBoundaryConditions();
    }
    else if (pimple_.dict().found("maxTInterfaceDdt"))
    {
        scalar maxDTdt
        (
            pimple_.dict().get<scalar>
            (
                "maxTInterfaceDdt"
            )
        );
        
        scalar dt(mesh_.time().deltaTValue());
        iT_ = saturation_->Tsat(p_);
        
        forAll(iT_, i)
        {
            const scalar& iT0(iT_.oldTime()[i]);
            scalar& iT(iT_[i]);
            scalar dTdt(min(max((iT-iT0)/dt, -maxDTdt), maxDTdt));
            iT = iT0 + dTdt*dt;
        }
    }
    else
    {
        iT_ = saturation_->Tsat(p_);
    }

    iT_.relax();
}


// ************************************************************************* //
