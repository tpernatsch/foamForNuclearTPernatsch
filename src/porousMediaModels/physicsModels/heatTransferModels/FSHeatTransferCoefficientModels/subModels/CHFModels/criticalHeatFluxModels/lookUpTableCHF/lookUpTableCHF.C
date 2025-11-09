/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2506                                                  |
|    Copyright 2011-2016 OpenFOAM Foundation, 2017-2025 OpenCFD Ltd.          |
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

#include "FSPair.H"
#include "FFPair.H"
#include "lookUpTableCHF.H"
#include "addToRunTimeSelectionTable.H"
#include "interpolation2DTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace criticalHeatFluxModels
{
    defineTypeNameAndDebug(lookUpTableCHF, 0);
    addToRunTimeSelectionTable
    (
        CHFModel,
        lookUpTableCHF,
        criticalHeatFluxModels
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::criticalHeatFluxModels::lookUpTableCHF::lookUpTableCHF
(
    const FSPair& pair,
    const dictionary& dict,
    const objectRegistry& objReg
)
:
    CHFModel
    (
        pair,
        dict,
        objReg
    ),
    mass_flow_quality_(pair.fluidRef().flowQuality()),
    p_(pair.mesh().lookupObject<volScalarField>("p")),
    rhoL_(pair.mesh().lookupObject<volScalarField>("thermo:rho.liquid")),
    rhoV_(pair.mesh().lookupObject<volScalarField>("thermo:rho.vapour")),
    uL_(pair.fluidRef().magU()),
    uV_(pair.mesh().lookupObject<volScalarField>("magU.vapour")),
    normalizedL_(pair.fluidRef().normalized()),
    FFPairPtr_(nullptr),
    Dh_(pair.structureRef().Dh()),
    PitchToDiameter_(dict.get<scalar>("PitchToDiameter")),
    pressureValues_(dict.lookup("pressureValues")),
    massFlowRateValues_(dict.lookup("massFlowRateValues")),
    qualityValues_(dict.lookup("qualityValues")),
    data_(PtrList<FieldField<Field, scalar>>(dict.lookup("data"), PtrListScalarFieldFieldINewGF())),
    pMethod(InterpolateTableBaseGF::interpolationMethodNames_[
        dict.lookupOrDefault<word>("pressureInterpolationMethod", "linear")
    ]),
    gMethod_(InterpolateTableBaseGF::interpolationMethodNames_[
        dict.lookupOrDefault<word>("massFlowRateInterpolationMethod", "linear")
    ]),
    xeMethod_(InterpolateTableBaseGF::interpolationMethodNames_[
        dict.lookupOrDefault<word>("qualityInterpolationMethod", "linear")
    ]),
    pTable_(pressureValues_, data_, pMethod)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::criticalHeatFluxModels::lookUpTableCHF::value
(
    const label& celli
) const
{

    // - Pointers breaking encapsulation - needed
    if (FFPairPtr_ == nullptr)
    {
        HashTable<const FFPair*> FFPairs(pair_.mesh().lookupClass<FFPair>());
        FFPairPtr_ = FFPairs[FFPairs.toc()[0]];
    }

    // ---------------------------------- //
    // -- Refs --//
    const scalar& xi(mass_flow_quality_[celli]);
    const scalar& pi(p_[celli]);
    const scalar& rhoLi(rhoL_[celli]);
    const scalar& rhoVi(rhoV_[celli]);
    const scalar& uLi(uL_[celli]);
    const scalar& uVi(uV_[celli]);
    const scalar& aLi(normalizedL_[celli]);
    const scalar& Dhi(Dh_[celli]);

    // Note that in this case, the fluid1 and fluid2 must be liquid and vapor, respectively.
    const fluid& fluid1(FFPairPtr_->fluid1());
    const fluid& fluid2(FFPairPtr_->fluid2());
    // Calculate the saturated variables field.
    const volScalarField& Tsat(FFPairPtr_->iT());
    volScalarField hLsat(fluid1.thermo().he(p_, Tsat));
    // Calculate the latent heat in this cell
    const scalar Li(mag(FFPairPtr_->L()[celli]));
    // Calculate the liquid and vapor saturated enthalpy in this cell.
    scalar hLsati(hLsat[celli]);
    scalar hVsati(hLsat[celli] + Li);
    // Calculate the mass flow rate in this cell
    scalar massFlowi(rhoLi*uLi*aLi+rhoVi*uVi*(1-aLi));
    // Calculate the liquid and vapor enthalpy in this cell.
    scalar hLi(fluid1.thermo().he()[celli]);
    scalar hVi(fluid2.thermo().he()[celli]);
    // Calculate the mixture enthalpy in this cell.
    // Since the enthalpy of vapour phase is always similar to saturated properties, we assume the enthalpy of vapour phase equal to saturated property.
//    scalar hMixi(xi*hLi+(1.0-xi)*hVi);
    scalar hMixi(xi*hLi+(1.0-xi)*hVsati);
    // Calculate the equilibrium quality in this cell
    scalar EquilibriumQualityi((hMixi-hLsati)/Li);

    //convert pressure to KPa, which is used in CHF look-up table
    scalar pk(pi/1000.0);

    // Interpolate 3D table with pressure
    FieldField<Field, scalar> pData(pTable_(pk));

    // Create 2D table (x is massFlowRate, y is quality)
    scalarFieldInterpolateTableGF gTable(massFlowRateValues_, pData, gMethod_);

    // Interpolate 2D table with mass flow rate
    scalarField gData(gTable(massFlowi));

    // Create 1D table (main coordinate is quality)
    scalarInterpolateTableGF xeTable(qualityValues_, gData, xeMethod_);

    // Interpolate 1D table and get final value of critical heat flux.
    // Then convert kW/m2 ----> W/m2
    scalar QCHF(xeTable(EquilibriumQualityi)*1e3);

    // Implement the correction
    // K1 for large diameter
    // K2 for rod bundle geometries
    if (PitchToDiameter_ != 1.0)
    {
      scalar K2((2.0*PitchToDiameter_ - 1.5)*exp(-pow(abs(EquilibriumQualityi),3.0)/2.0));
      QCHF = QCHF * K2;
    }
    else
    {
      scalar K1(max(0.6,sqrt(0.008/Dhi)));
      QCHF = QCHF * K1;
    }

    return QCHF;
}