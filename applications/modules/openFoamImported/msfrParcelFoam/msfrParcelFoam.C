/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2512                                                  |
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

#include "msfrParcelFoam.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(msfrParcelFoam, 0);
    addToRunTimeSelectionTable
    (
        solver,
        msfrParcelFoam,
        dynamicFvMesh
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::msfrParcelFoam::msfrParcelFoam
(
    dynamicFvMesh& mesh
)
:
    solver(mesh),
    pimple_(mesh_),
    npimple_(mesh_, "NPIMPLE"),
    U_
    (
        IOobject
        (
            "U",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    phi_
    (
        IOobject
        (
            "phi",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::flux(U_)
    ),
    Uf_
    (
        new surfaceVectorField
        (
            IOobject
            (
                "Uf",
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fvc::interpolate(U_)
        )
    ),
    laminarTransport_(U_, phi_),
    beta_
    (
        "beta",
        dimless/dimTemperature,
        laminarTransport_
    ),
    CpRef_
    (
        "CpRef",
        dimEnergy/dimMass/dimTemperature,
        laminarTransport_
    ),
    TRef_
    (
        "TRef",
        dimTemperature,
        laminarTransport_
    ),
    Pr_
    (
        "Pr",
        dimless,
        laminarTransport_
    ),
    Prt_
    (
        "Prt",
        dimless,
        laminarTransport_
    ),
    Sc_
    (
        "Sc",
        dimless,
        laminarTransport_
    ),
    Sct_
    (
        "Sct",
        dimless,
        laminarTransport_
    ),
    T_
    (
        IOobject
        (
            "T",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    p_rgh_
    (
        IOobject
        (
            "p_rgh",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    vort_
    (
        IOobject
        (
            "vort",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("zero", dimensionSet(0, 0, -1, 0, 0, 0, 0), vector::zero)
    ),
    turbulence_
    (
        incompressible::turbulenceModel::New(U_, phi_, laminarTransport_)
    ),
    rhoInfValue_
    (
        "rhoInf",
        dimDensity,
        laminarTransport_
    ),
    rhoInf_
    (
        IOobject
        (
            "rho",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        rhoInfValue_
    ),
    muc_
    (
        IOobject
        (
            "muc",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rhoInf_*laminarTransport_.nu()
    ),
    rhok_
    (
        IOobject
        (
            "rhok",
            mesh_.time().timeName(),
            mesh_
        ),
        1.0 - beta_*(T_ - TRef_)
    ),
    alphat_
    (
        IOobject
        (
            "alphat",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    g_
    (
        IOobject
        (
            "g",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    hRef_
    (
        IOobject
        (
            "hRef",
            mesh_.time().constant(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        dimensionedScalar("hRef", dimLength, 0)
    ),
    gh_
    (
        "gh",
        (g_ & mesh_.C()) + mag(g_)*hRef_
    ),
    ghf_
    (
        "ghf",
        (g_ & mesh_.Cf()) + mag(g_)*hRef_
    ),
    p_
    (
        IOobject
        (
            "p",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p_rgh_ + rhok_*gh_
    ),
    pRefCell_(0),
    pRefValue_(0.0),
    Qdot_
    (
        IOobject
        (
            "Qdot",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero)
    ),
    MRF_(mesh_),
    radiation_(radiation::radiationModel::New(T_)),
    rhoCpRef_
    (
        "rhoCpRef",
        rhoInfValue_*CpRef_
    ),
    fvOptions_(fv::options::New(mesh_)),
    FP_
    (
        "FPCloud",
        rhoInf_,
        U_,
        muc_,
        vort_,
        T_,
        g_
    ),
    tsurfaceFilm_(regionModels::surfaceFilmModel::New(mesh_, g_, "surfaceFilm")),
    surfaceFilm_(tsurfaceFilm_()),
    solvePrimaryRegion_(pimple_.dict().getOrDefault("solvePrimaryRegion", true)),
    cumulativeContErr_(0),

    // --- Neutronics / decay heat / fission product transport ---

    nuclearPropertiesDict_
    (
        IOobject
        (
            "nuclearProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    nEnergyGroups_
    (
        readInt(nuclearPropertiesDict_.subDict("neutronTransport").lookup("energyGroups"))
    ),
    useDoppler_(nuclearPropertiesDict_.subDict("neutronTransport").lookup("useDoppler")),
    DopplerControl_("DopplerControl", dimless, useDoppler_ ? 1 : 0),
    rhoRefXS_("rhoRefXS", dimDensity, nuclearPropertiesDict_.subDict("neutronTransport")),
    TRefXS_("TRefXS", dimTemperature, nuclearPropertiesDict_.subDict("neutronTransport")),
    D_ref_(nEnergyGroups_),
    alphaD_ref_(nEnergyGroups_),
    Sa_ref_(nEnergyGroups_),
    alphaSa_ref_(nEnergyGroups_),
    Sf_ref_(nEnergyGroups_),
    alphaSf_ref_(nEnergyGroups_),
    Ss_ref_(nEnergyGroups_*nEnergyGroups_),
    alphaSs_ref_(nEnergyGroups_*nEnergyGroups_),
    Nu_(nEnergyGroups_),
    Ef_(nEnergyGroups_),
    invVel_(nEnergyGroups_),
    chiPrompt_(nEnergyGroups_),
    chiDelayed_(nEnergyGroups_),
    nPrecGroups_
    (
        readInt(nuclearPropertiesDict_.subDict("delayedNeutronPrecursors").lookup("groups"))
    ),
    precLambda_(nPrecGroups_),
    precBeta_(nPrecGroups_),
    precBetaTot_("precBetaTot", dimless, 0),
    nDecGroups_
    (
        readInt(nuclearPropertiesDict_.subDict("decayHeatPrecursors").lookup("groups"))
    ),
    decLambda_(nDecGroups_),
    decBeta_(nDecGroups_),
    decBetaTot_("decBetaTot", dimless, 0),
    fpTransport_(nuclearPropertiesDict_.found("fissionProducts")),
    nFPSpecies_(0),
    flux_(nEnergyGroups_),
    fluxTot_
    (
        IOobject
        (
            "fluxTot",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimless/dimArea/dimTime, 0)
    ),
    prec_(nPrecGroups_),
    dec_(nDecGroups_),
    D_(nEnergyGroups_),
    Sa_(nEnergyGroups_),
    Sf_(nEnergyGroups_),
    Ss_(nEnergyGroups_*nEnergyGroups_),
    M_(nEnergyGroups_),
    unitField_
    (
        IOobject
        (
            "unitField",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimless, 1)
    ),
    zeroField_("zeroField", 0*unitField_),
    Rf_
    (
        IOobject
        (
            "Rf",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimless/dimVolume/dimTime, 0)
    ),
    q_
    (
        IOobject
        (
            "q",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimPower/dimVolume, 0)
    ),
    qPrompt_
    (
        IOobject
        (
            "qPrompt",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimPower/dimVolume, 0)
    ),
    qDecay_
    (
        IOobject
        (
            "qDecay",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimPower/dimVolume, 0)
    ),
    Q_
    (
        IOobject
        (
            "Q",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimPower, 0)
    ),
    Qprompt_
    (
        IOobject
        (
            "Qprompt",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimPower, 0)
    ),
    Qdecay_
    (
        IOobject
        (
            "Qdecay",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimPower, 0)
    ),
    NuRf_
    (
        IOobject
        (
            "NuRf",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimless/dimVolume/dimTime, 0)
    ),
    EfRf_
    (
        IOobject
        (
            "EfRf",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimEnergy/dimVolume/dimTime, 0)
    ),
    neutronSource_(nEnergyGroups_),
    rhoNorm_
    (
        IOobject
        (
            "rhoNorm",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimless, 0)
    ),
    logT_("logT", zeroField_),
    tempCoupling_(mesh_.time().controlDict().lookup("tempCoupling")),
    criticality_(mesh_.time().controlDict().lookupOrDefault<word>("criticality", "yes")),
    KeffInit_
    (
        "KeffInit",
        dimless,
        mesh_.time().controlDict().lookupOrDefault<scalar>("Keff", 1)
    ),
    Qnominal_("nominalPower", dimPower, mesh_.time().controlDict()),
    Keff_
    (
        IOobject
        (
            "Keff",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        KeffInit_
    ),
    reactivity_
    (
        IOobject
        (
            "reactivity",
            mesh_.time().timeName(),
            mesh_
        ),
        (Keff_ - 1)/Keff_
    )
{
    turbulence_->validate();

    setRefCell
    (
        p_,
        p_rgh_,
        pimple_.dict(),
        pRefCell_,
        pRefValue_
    );

    if (p_rgh_.needReference())
    {
        p_ += dimensionedScalar
        (
            "p",
            p_.dimensions(),
            pRefValue_ - getRefCellValue(p_, pRefCell_)
        );
    }

    mesh_.setFluxRequired(p_rgh_.name());

    Info<< "    Nominal reactor power  : " << Qnominal_.value()/1E+06 << " MW" << endl;
    Info<< "    Keff                   : " << max(Keff_).value() << endl;
    Info<< "    Reactivity             : " << max(reactivity_).value()*1e5 << " pcm" << endl;

    #include "include/create/readNuclearProperties.H"

    if (fpTransport_)
    {
        #include "include/create/readFissionProductsProperties.H"
    }

    #include "include/create/createNuclearFields.H"

    if (fpTransport_)
    {
        #include "include/create/createFissionProductsFields.H"
    }

    #include "include/equations/updateCrossSections.H"
    #include "include/equations/updateFissionRate.H"
    #include "include/equations/updatePowerSource.H"
    #include "include/equations/updateNeutronSource.H"
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::msfrParcelFoam::correctPhysics()
{
    FP_.storeGlobalPositions();
    FP_.evolve();

    surfaceFilm_.evolve();

    if (solvePrimaryRegion_)
    {
        while (pimple_.loop())
        {
            #include "include/equations/UEqn.H"

            vort_ = fvc::curl(U_);

            #include "include/equations/TEqn.H"

            while (pimple_.correct())
            {
                #include "include/equations/pEqn.H"
            }

            if (pimple_.turbCorr())
            {
                laminarTransport_.correct();
                turbulence_->correct();
            }

            #include "include/equations/updateCrossSections.H"

            while (npimple_.loop())
            {
                #include "include/equations/fluxEqns.H"
                #include "include/equations/precEqns.H"
                #include "include/equations/decEqns.H"

                #include "include/equations/updateNeutronSource.H"
                #include "include/equations/updateFissionRate.H"
            }

            #include "include/equations/updatePowerSource.H"

            #include "include/equations/fpEqns.H"
        }
    }

    Info<< "Q (MW)      : " << max(Q_).value()/1E+06 << endl;
}


Foam::scalar Foam::solvers::msfrParcelFoam::maxDeltaT()
{
    scalar newDeltaT = mesh_.time().controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);

    scalar maxCo = mesh_.time().controlDict().getOrDefault<scalar>("maxCo", 1);

    if (mesh_.time().value() > mesh_.time().controlDict().get<scalar>("deltaT"))
    {
        scalar CoNum = 0.0;
        scalar meanCoNum = 0.0;

        if (mesh_.nInternalFaces())
        {
            scalarField sumPhi(fvc::surfaceSum(mag(phi_))().primitiveField());

            CoNum = 0.5*gMax(sumPhi/mesh_.V().field())*mesh_.time().deltaTValue();

            meanCoNum =
                0.5*(gSum(sumPhi)/gSum(mesh_.V().field()))*mesh_.time().deltaTValue();
        }

        Info<< "Courant Number mean: " << meanCoNum
            << " max: " << CoNum << endl;

        scalar maxDeltaTFact = maxCo/(CoNum + SMALL);
        scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

        newDeltaT = min(deltaTFact*mesh_.time().deltaTValue(), newDeltaT);
    }

    return newDeltaT;
}


// ************************************************************************* //
