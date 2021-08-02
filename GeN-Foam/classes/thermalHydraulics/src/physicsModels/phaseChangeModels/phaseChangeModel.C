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
#include "FFPair.H"
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
    FFPair& pair,
    const dictionary& dict
)
:
    IOdictionary
    (
        IOobject
        (
            IOobject::groupName("phaseChangeModel", typeName),
            pair.mesh().time().constant(),
            pair.mesh()
        ),
        dict
    ),
    mesh_(pair.mesh()),
    pimple_(pair.pimple()),
    pair_(pair),
    fluid1_(pair.fluid1()),
    fluid2_(pair.fluid2()),
    liquid_
    (
        (pair.fluid1().isLiquid()) ? 
        pair.fluid1() : pair.fluid2()
    ),
    vapour_
    (
        (pair.fluid1().isGas()) ? 
        pair.fluid1() : pair.fluid2()
    ),
    p_(mesh_.lookupObject<volScalarField>("p")),
    T1_(fluid1_.thermo().T()),
    T2_(fluid2_.thermo().T()),
    htc1_(pair.htc(fluid1_.name())),
    htc2_(pair.htc(fluid2_.name())),
    iT_(pair.iT()),
    iA_(pair.iA()),
    L_
    (
        IOobject
        (
            IOobject::groupName("L", pair.name()),
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimEnergy/dimMass, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    dmdt_
    (
        IOobject
        (
            IOobject::groupName("dmdt", pair.name()),
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimDensity/dimTime, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    latentHeatPtr_
    (
        latentHeatModel::New
        (
            *this,
            this->subDict("latentHeatModel"),
            mesh_
        )
    ),
    saturationPtr_
    (
        saturationModel::New
        (
            *this,
            this->subDict("saturationModel"),
            mesh_
        )
    ),
    residualIA_
    (
        dict.getOrDefault<scalar>("residualInterfacialArea", 1e-3)
    ),
    residualIACells_(0),
    dmLostToLimiter_(0)
{
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
        saturationPtr_->correctField(iT_, "TSat");
    }

    //-
    latentHeatPtr_->correctField(L_);

    //- Read residualIACells from cellZones
    wordList residualIARegions
    (
        this->lookupOrDefault<wordList>
        (
            "residualInterfacialAreaRegions", 
            wordList()
        )
    );
    forAll(residualIARegions, i)
    {
        const labelList& regionCells
        (
            mesh_.cellZones()[residualIARegions[i]]
        );
        forAll(regionCells, j)
        {
            residualIACells_.append(regionCells[j]);
        }
    }
}

// * * * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * //

void Foam::phaseChangeModel::limitInterfacialArea()
{
    if (residualIA_ != 0.0)
    {
        if (residualIACells_.size() > 0)
        {
            forAll(residualIACells_, i)
            {
                const label& celli(residualIACells_[i]);
                scalar& iAi(iA_[celli]);
                iAi = max(iAi, residualIA_);
            }
        }
        else
        {
            forAll(mesh_.cells(), celli)
            {
                scalar& iAi(iA_[celli]);
                iAi = max(iAi, residualIA_);
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

        /*if (pimple_.finalIter())
        {
            Info<< "Cumulative dm lost to limiter = " 
                << (dmLostToLimiter_) 
                << " kg/m3" << endl;
        }*/
    }

    dmdt_.correctBoundaryConditions();
}

void Foam::phaseChangeModel::correctInterfacialTemperature()
{
    if (mesh_.relaxField("T.interface"))
        iT_.storePrevIter();
    if (pimple_.dict().found("maxTInterfaceDdt"))
    {
        scalar maxDTdt
        (
            pimple_.dict().get<scalar>
            (
                "maxTInterfaceDdt"
            )
        );
        
        scalar dt(mesh_.time().deltaTValue());
        saturationPtr_->correctField(iT_, "TSat");
        
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
        saturationPtr_->correctField(iT_, "TSat");
    }

    iT_.relax();
}

// ************************************************************************* //
