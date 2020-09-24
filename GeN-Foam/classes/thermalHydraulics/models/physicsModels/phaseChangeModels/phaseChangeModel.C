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
    residualIACells_(0),
    dmLostToLimiter_(0)
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

    /*
    All right, what's the deal with this apparent BS I am doing here? The thing
    is, over countless simulations (of the KNS-37 case in 2D) I noticed that
    ultimate root cause of all crashed simulations consisted in spikes in the
    vapour temperature. While a (severe) under-relaxation of the vapour
    temperature saves a considerable amount of cases, occasionaly crashes would
    still occur far after the bulk of the time dynamics of the transient was
    resolved. Always this was observed to be due to/result in vapour 
    temperature spikes. Now, I wrote due to/result in as it is difficult to
    figure out which is the first parameter to fuck away given the tight 
    coupling between all of the physics. However, ultimately, everything seems
    to start to mass transfer spikes in response to continity error spikes, 
    which seem to come first after all. Thus, here is the idea: limit
    dmdt in those cells in which the (normalized) continuity error of a phase
    is above a maximum or below a minimum allowed value. Given that
    dmdt is positive in those cell where phase 1 is disappearing and becoming
    phase 2, a maximum dmdt value can be fixed in those cells so to avoid
    the enitre phase fraction to disappear in a single time step. On the flip
    side, dmdt is negative in those cells where phase 2 is disappearing and
    becoming phase 1, so a minimum value of dmdt can be fixed based on the same
    logic. So, given that the continuity equations are:

    ddt(alpha1*rho1) + div(alphaRhoPhi1) = -dmdt
    ddt(alpha2*rho2) + div(alphaRhoPhi2) = dmdt

    assuming the entire phase fraction to disappear in a single time step
    allows us to compute the minum allowable ddt(alpha_i) as 
    -alpha_i*rho_i/deltaT. The maximum value of dmdt is obtained by the 
    continuity eq of phase 1 and the minum from the eq of phase 2.
    In particular, the maximum is fixed only in those cells where the 
    normalized continuity error of phase 1 is above a certain threshold
    selected by the used. Specularly, the minimum is fixed only in those
    cells where the normalized continuity error of phase 2 is above a certain
    threshold. This is user-selectable, values in the 0.1-1 range reasonable
    (i.e. cells in which the continuity errors are from 10% to 100% of the 
    total fluid mass in the cell... yes, they do unfortunately occur during
    violent transients if the time steps are relatively large)
    */
    if 
    (
        pimple_.dict().lookupOrDefault<bool>
        (
            "adaptiveMassTransferLimiter", 
            false
        )
    )
    {
        scalar SF0
        (
            pimple_.dict().lookupOrDefault<scalar>
            (
                "adaptiveMassTransferSafetyFactor", 1.0
            )
        );
        scalar maxFrac
        (
            pimple_.dict().lookupOrDefault<scalar>
            (
                "adaptiveMassTransferMaxAdjustmentRatio", 2.0
            )
        );
        
        const scalarField& V(mesh_.V());
        scalar totV(0);
        forAll(V, i)
        {
            totV += V[i];
        }
        scalar dt(mesh_.time().deltaT().value());

        scalarField div1(fvc::div(fluid1_.alphaRhoPhi())().primitiveField());
        scalarField div2(fvc::div(fluid2_.alphaRhoPhi())().primitiveField());
        scalarField contErrNorm1
        (
            mag
            (
                fluid1_.contErr().primitiveField()/
                fluid1_.rho().primitiveField()
            )
        );
        scalarField contErrNorm2
        (
            mag
            (
                fluid2_.contErr().primitiveField()/
                fluid2_.rho().primitiveField()
            )
        );
        
        scalar minContErr(1e-6);
        
        scalar posDMdt(0.0);
        scalar posLimitedDMdt(0.0);
        scalar negDMdt(0.0);
        scalar negLimitedDMdt(0.0);
        
        scalar dmdtAvgBeforeLim(0.0);
        if (pimple_.finalIter())
        {
            dmdtAvgBeforeLim = dmdt_.weightedAverage(mesh_.V()).value();
        }

        forAll(dmdt_, i)
        {
            scalar& dmdt(dmdt_[i]);
            if (dmdt > 0.0)
            {
                const scalar& Vi(V[i]);
                scalar SF(SF0/max(contErrNorm1[i], minContErr));
                scalar maxDmdt
                (
                    max(SF*fluid1_[i]*fluid1_.rho()[i]/dt-div1[i], 0.0)
                );
                posDMdt += dmdt*Vi;
                dmdt = min(dmdt, maxDmdt);
                posLimitedDMdt += dmdt*Vi;
            }
            else if (dmdt < 0.0)
            {
                const scalar& Vi(V[i]);
                scalar SF(SF0/max(contErrNorm2[i], minContErr));
                scalar minDmdt
                (
                    min(-SF*fluid2_[i]*fluid2_.rho()[i]/dt+div2[i], 0.0)
                );
                negDMdt += dmdt*Vi;
                dmdt = max(dmdt, minDmdt);
                negLimitedDMdt += dmdt*Vi;
            }
        }

        //- Re-absorb lost dm in total dmdt integral as if it were to be
        //  compensated in a single time step. I don;t know if it makes sense
        //  to reintroduce what was lost precisely because I assume the dmdt
        //  in high continuity error cells should not be trusted... oh, well!
        //  Defaults to false anyway
        if 
        (
            pimple_.dict().lookupOrDefault<bool>
            (
                "adaptiveMassTransferReinsertion", 
                false
            )
        )
        {
            scalar dMdtLoss(dmLostToLimiter_*totV/dt);
            if (dMdtLoss > 0.0)
                negDMdt -= dMdtLoss;
            else
                posDMdt -= dMdtLoss;
        }

        scalar posFrac
        (
            max(min(posDMdt/max(posLimitedDMdt, 1e-69), maxFrac), 1.0)
        );
        scalar negFrac
        (
            max(min(negDMdt/min(negLimitedDMdt, -1e-69), maxFrac), 1.0)
        );  
        forAll(dmdt_, i)
        {
            scalar& dmdt(dmdt_[i]);
            if (dmdt > 0.0)
            {
                dmdt *= posFrac;
            }
            else if (dmdt < 0.0)
            {
                dmdt *= negFrac;
            }
        }
           
        if (pimple_.finalIter())
        {
            scalar dmdtAvgAfterLim(dmdt_.weightedAverage(mesh_.V()).value());
            dmLostToLimiter_ += (dmdtAvgAfterLim - dmdtAvgBeforeLim)*dt;
            Info<< "Cumulative dm lost to limiter = " 
                << (dmLostToLimiter_) 
                << " kg/m3" << endl;
        }
    }

    dmdt_.correctBoundaryConditions();
}


void Foam::phaseChangeModel::correctInterfacialT()
{
    iT_.storePrevIter();

    if (pimple_.dict().found("maxTInterfaceIterRelChange"))
    {
        scalar maxRelChange
        (
            pimple_.dict().get<scalar>
            (
                "maxTInterfaceIterRelChange"
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
            const scalar& iT0(iT_.prevIter()[i]);
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
