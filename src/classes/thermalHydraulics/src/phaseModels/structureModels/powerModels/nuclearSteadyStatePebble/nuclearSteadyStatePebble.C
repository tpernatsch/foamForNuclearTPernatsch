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

#include "nuclearSteadyStatePebble.H"
#include "structure.H"
#include "addToRunTimeSelectionTable.H"
#include "SquareMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace powerModels
{
    defineTypeNameAndDebug(nuclearSteadyStatePebble, 0);
    addToRunTimeSelectionTable
    (
        powerModel,
        nuclearSteadyStatePebble,
        powerModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::powerModels::nuclearSteadyStatePebble::nuclearSteadyStatePebble
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
    fractionOfPowerFromNeutronics_(0),
    TpS_
    (
        IOobject
        (
            "TpS."+typeName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    Tmav_
    (
        IOobject
        (
            "Tmav."+typeName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    Tmout_
    (
        IOobject
        (
            "Tmout."+typeName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    TfS_
    (
        IOobject
        (
            "TfS."+typeName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    Tfav_
    (
        IOobject
        (
            "Tfav."+typeName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    Tfmax_
    (
        IOobject
        (
            "Tfmax."+typeName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    Tfavmax_(0),
    Tfavmin_(1e69),
    Tfmaxmax_(0),
    Tfmaxmin_(1e69),
    Tmoutmax_(0),
    Tmoutmin_(1e69),
    Tmavmax_(0),
    Tmavmin_(1e69),
    TpSmax_(0),
    TpSmin_(1e69),
    TfSmax_(0),
    TfSmin_(1e69),
    rcore_(0),     // Pebble inner graphite core
    rmatrix_(0),   // Pebble fuel layer with matrix
    rshell_(0),    // Pebble outer graphite shell
    rfuel_(0),     // TRISO fuel kernel
    rbuffer_(0),   // TRISO buffer layer
    rPyCin_(0),    // TRISO inner PyC layer
    rSiC_(0),      // TRISO SiC layer
    rPyCout_(0),   // TRISO outer PyC layer
    kfuel_(0),     // TRISO fuel kernel
    kbuffer_(0),   // TRISO buffer layer
    kPyC_(0),    // TRISO PyC layers
    kSiC_(0),      // TRISO SiC layer
    kgraphite_(0), // Pebble matrix graphite
    keffMatrix_
    (
        IOobject
        (
            "keffmatrix."+typeName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimPower/dimLength/dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    rhofuel_(
        IOobject
        (
            "rhofuel."+typeName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimMass/dimVolume, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    rhobuffer_(
        IOobject
        (
            "rhobuffer."+typeName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimMass/dimVolume, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    rhoPyC_(
        IOobject
        (
            "rhoPyC."+typeName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimMass/dimVolume, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    rhoSiC_(
        IOobject
        (
            "rhoSiC."+typeName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimMass/dimVolume, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    rhographite_(
        IOobject
        (
            "rhographite."+typeName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimMass/dimVolume, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    rhocoolant_(
        IOobject
        (
            "rhocoolant."+typeName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimMass/dimVolume, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    surfaceHeatFlux_
    (
        IOobject
        (
            "surfaceHeatFlux."+typeName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimPower/dimLength/dimLength, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    Ntriso_(0),   // Number of TRISO per pebble
    cellToRegion_(mesh_.cells().size(), 0),
    regionIndexToRegionName_(0),
    pi_(constant::mathematical::pi)
{
    structure_.setRegionField(*this, structureRef.powerDensityNeutronics(), "powerDensityNeutronics");
    const scalarField& V(mesh_.V());
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
        scalar fractionOfPowerFromNeutronics(dict.lookupOrDefault<scalar>("fractionOfPowerFromNeutronics",1.0));
        scalar TpS_guess(dict.lookupOrDefault<scalar>("pebbleSurfaceTemperatureGuess", 900));
        scalar rcore(dict.get<scalar>("pebbleCoreRadius"));
        scalar rmatrix(dict.get<scalar>("pebbleMatrixRadius"));
        scalar rshell(dict.get<scalar>("pebbleShellRadius"));
        scalar rfuel(dict.get<scalar>("trisoFuelRadius"));
        scalar rbuffer(dict.get<scalar>("trisoBufferRadius"));
        scalar rPyCin(dict.get<scalar>("trisoInnerPyCRadius"));
        scalar rSiC(dict.get<scalar>("trisoSiCRadius"));
        scalar rPyCout(dict.get<scalar>("trisoOuterPyCRadius"));

        label  Ntriso(dict.get<label>("nTRISO"));

        Polynomial<8> kfuelCoeffs(dict.get<Polynomial<8>>("trisoFuelKCoeffs"));
        Polynomial<8> kbufferCoeffs(dict.get<Polynomial<8>>("trisoBufferKCoeffs"));
        Polynomial<8> kSiCCoeffs(dict.get<Polynomial<8>>("trisoSiCKCoeffs"));
        Polynomial<8> kPyCCoeffs(dict.get<Polynomial<8>>("trisoPyCKCoeffs"));
        Polynomial<8> kgraphiteCoeffs(dict.get<Polynomial<8>>("pebbleGraphiteKCoeffs"));

        Polynomial<8> rhofuelCoeffs(dict.get<Polynomial<8>>("trisoFuelDensityCoeffs"));
        Polynomial<8> rhobufferCoeffs(dict.get<Polynomial<8>>("trisoBufferDensityCoeffs"));
        Polynomial<8> rhoSiCCoeffs(dict.get<Polynomial<8>>("trisoSiCDensityCoeffs"));
        Polynomial<8> rhoPyCCoeffs(dict.get<Polynomial<8>>("trisoPyCDensityCoeffs"));
        Polynomial<8> rhographiteCoeffs(dict.get<Polynomial<8>>("pebbleGraphiteDensityCoeffs"));
        Polynomial<8> rhocoolantCoeffs(dict.get<Polynomial<8>>("coolantDensityCoeffs"));

        //- Fill in lists for this region

        fractionOfPowerFromNeutronics_.append(fractionOfPowerFromNeutronics);
        TpS_guess_.append(TpS_guess);

        rcore_.append(rcore);
        rmatrix_.append(rmatrix);
        rshell_.append(rshell);
        rfuel_.append(rfuel);
        rbuffer_.append(rbuffer);
        rPyCin_.append(rPyCin);
        rSiC_.append(rSiC);
        rPyCout_.append(rPyCout);

        Ntriso_.append(Ntriso);

        kfuelCoeffs_.append(kfuelCoeffs);
        kbufferCoeffs_.append(kbufferCoeffs);
        kSiCCoeffs_.append(kSiCCoeffs);
        kPyCCoeffs_.append(kPyCCoeffs);
        kgraphiteCoeffs_.append(kgraphiteCoeffs);

        rhofuelCoeffs_.append(rhofuelCoeffs);
        rhobufferCoeffs_.append(rhobufferCoeffs);
        rhoSiCCoeffs_.append(rhoSiCCoeffs);
        rhoPyCCoeffs_.append(rhoPyCCoeffs);
        rhographiteCoeffs_.append(rhographiteCoeffs);
        rhocoolantCoeffs_.append(rhocoolantCoeffs);
    }

    forAll(this->cellList_, i)
    {
        label celli(this->cellList_[i]);

        label regioni(cellToRegion_[celli]);
        scalar rcore(rcore_[regioni]);
        scalar rmatrix(rmatrix_[regioni]);
        scalar rshell(rshell_[regioni]);
        scalar rfuel(rfuel_[regioni]);
        scalar rbuffer(rbuffer_[regioni]);
        scalar rPyCin(rPyCin_[regioni]);
        scalar rSiC(rSiC_[regioni]);
        scalar rPyCout(rPyCout_[regioni]);

        label  Ntriso(Ntriso_[regioni]);

        scalar TpS_guess(TpS_guess_[regioni]);

        scalar kfuel(kfuelCoeffs_[regioni].value(TpS_guess));
        scalar kbuffer(kbufferCoeffs_[regioni].value(TpS_guess));
        scalar kSiC(kSiCCoeffs_[regioni].value(TpS_guess));
        scalar kPyC(kPyCCoeffs_[regioni].value(TpS_guess));
        scalar kgraphite(kgraphiteCoeffs_[regioni].value(TpS_guess));

        scalar rhofuel(rhofuelCoeffs_[regioni].value(TpS_guess));
        scalar rhobuffer(rhobufferCoeffs_[regioni].value(TpS_guess));
        scalar rhoSiC(rhoSiCCoeffs_[regioni].value(TpS_guess));
        scalar rhoPyC(rhoPyCCoeffs_[regioni].value(TpS_guess));
        scalar rhographite(rhographiteCoeffs_[regioni].value(TpS_guess));
        scalar rhocoolant(rhocoolantCoeffs_[regioni].value(TpS_guess));

        // Initial calculation of the temperatures
        scalar q = structure_.powerDensityNeutronics()[celli];
        const scalar& dV(V[celli]);
        scalar Q_cell = q*dV;
        scalar Q_pebble_tot =
            alpha_[celli]== 0 ?
            0
            :Q_cell*((4.0/3.0 * pi_ * pow(rshell, 3))/dV)*(1.0/alpha_[celli]);
        //scalar q_pebble = Q_pebble_tot / (4.0/3.0 * pi_ * pow(rshell, 3));
        scalar q_matrix = Q_pebble_tot / (4.0/3.0 * pi_ * (pow(rmatrix, 3) - pow(rcore, 3)));
        scalar Q_triso = Q_pebble_tot / Ntriso;
        scalar q_fuel = Q_triso / (4.0/3.0 * pi_ * pow(rfuel, 3));
        surfaceHeatFlux_[celli] = Q_pebble_tot / (4.0 * pi_ * pow(rshell, 2));

        // Calculation of effective conductivity in matrix (Maxwell model)
        scalar ktriso = (1/rfuel-1/rPyCout)/((1/kbuffer)*(1/rfuel-1/rbuffer)+(1/kPyC)*(1/rbuffer-1/rPyCin)+(1/kSiC)*(1/rPyCin-1/rSiC)+1/kPyC*(1/rSiC-1/rPyCout));
        scalar PF = (4.0/3.0*pi_*pow(rPyCout,3)) * static_cast<double >(Ntriso) / (4.0/3.0*pi_*(pow(rmatrix,3)-pow(rcore,3)));
        scalar kappa_coef = ktriso/kgraphite;
        scalar beta_coef  = (kappa_coef-1)/(kappa_coef+2);
        scalar keffmatrix = kgraphite*((1+2*beta_coef*PF)/(1-beta_coef*PF));
        //Info << PF << " " << ktriso << " " << kappa_coef << " " << beta_coef << " " << keffmatrix << " " << kfuel << " " << kbuffer << " " << kSiC << " " << kPyC << " " << kgraphite << endl;
        // Calculation of average matrix and fuel temperatures directly
        TpS_[celli]   = TpS_guess;
        Tmout_[celli] = TpS_[celli] + (Q_pebble_tot/(4*pi_*kgraphite))*(1/rmatrix-1/rshell);
        Tmav_[celli]  = Tmout_[celli]+(q_matrix/(6*keffmatrix))*(pow(rmatrix,2)-3*(pow(rmatrix,5)-pow(rcore,5))/(5*(pow(rmatrix,3)-pow(rcore,3))))+(pow(rcore,3)*q_matrix)/(3*keffmatrix)*(1/rmatrix-3*(pow(rmatrix,2)-pow(rcore,2))/(2*(pow(rmatrix,3)-pow(rcore,3))));
        TfS_[celli]   = Tmav_[celli] + Q_triso/(4*pi_) * (1/kbuffer*(1/rfuel-1/rbuffer) +  1/kPyC*(1/rbuffer-1/rPyCin) + 1/kSiC*(1/rPyCin-1/rSiC) + 1/kPyC*(1/rSiC-1/rPyCout));
        Tfav_[celli]  = TfS_[celli] + q_fuel*pow(rfuel,2) / (15*kfuel);
        Tfmax_[celli]  = TfS_[celli] + q_fuel*pow(rfuel,2) / (6*kfuel);
        keffMatrix_[celli] = keffmatrix;
        rhofuel_[celli] = rhofuel;
        rhobuffer_[celli] = rhobuffer;
        rhoSiC_[celli] = rhoSiC;
        rhoPyC_[celli] = rhoPyC;
        rhographite_[celli] = rhographite;
        rhocoolant_[celli] = rhocoolant;
    }

    //- Initialize in dict
    this->IOdictionary::set("Tfav", Tfav_);
    this->IOdictionary::set("Tfmax", Tfmax_);
    this->IOdictionary::set("Tmav", Tmav_);
    this->IOdictionary::set("TpS", TpS_);
    this->IOdictionary::set("Tmout", Tmout_);
    this->IOdictionary::set("TfS", TfS_);
    this->IOdictionary::set("keffmatrix", keffMatrix_);
    this->IOdictionary::set("surfaceHeatFlux", surfaceHeatFlux_);
    this->IOdictionary::set("rhofuel", rhofuel_);
    this->IOdictionary::set("rhobuffer", rhobuffer_);
    this->IOdictionary::set("rhoSiC", rhoSiC_);
    this->IOdictionary::set("rhoPyC", rhoPyC_);
    this->IOdictionary::set("rhographite", rhographite_);
    this->IOdictionary::set("rhocoolant", rhocoolant_);

    Tfav_.correctBoundaryConditions();
    Tfmax_.correctBoundaryConditions();
    Tmav_.correctBoundaryConditions();
    TpS_.correctBoundaryConditions();
    Tmout_.correctBoundaryConditions();
    TfS_.correctBoundaryConditions();
    rhofuel_.correctBoundaryConditions();
    rhobuffer_.correctBoundaryConditions();
    rhoSiC_.correctBoundaryConditions();
    rhoPyC_.correctBoundaryConditions();
    rhographite_.correctBoundaryConditions();
    rhocoolant_.correctBoundaryConditions();
    keffMatrix_.correctBoundaryConditions();
    surfaceHeatFlux_.correctBoundaryConditions();

/* The equivalent of this should be added
        // Update average fuel and clad temp used for coupling
        this->structureRef().TFuelAv()[celli] = Tfav_[celli];
        this->structureRef().TCladAv()[celli] = Tcav_[celli];
*/
    //- Finally, set up interfacial area
    this->setInterfacialArea();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::powerModels::nuclearSteadyStatePebble::~nuclearSteadyStatePebble()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::powerModels::nuclearSteadyStatePebble::setInterfacialArea()
{
    forAll(this->cellList_, i)
    {
        const label& celli(this->cellList_[i]);
        iA_[celli] = 3.0*alpha_[celli]/(rshell_[cellToRegion_[celli]]);
    }
    iA_.correctBoundaryConditions();
}

void Foam::powerModels::nuclearSteadyStatePebble::correct
(
    const volScalarField& HTSum,  // == SUM_j [htc_j*T_j*frac_j]
    const volScalarField& HSum    // == SUM_j [htc_j*frac_j]
)
{
    //- Reset min, max, average temperatures
    Tfavmax_ = 0.0;
    Tfavmin_ = 1e69;
    Tfmaxmax_ = 0.0;
    Tfmaxmin_ = 1e69;
    Tmoutmax_ = 0.0;
    Tmoutmin_ = 1e69;
    Tmavmax_ = 0.0;
    Tmavmin_ = 1e69;
    TpSmax_ = 0.0;
    TpSmin_ = 1e69;
    TfSmax_ = 0.0;
    TfSmin_ = 1e69;

    const scalarField& V(mesh_.V());
    scalar totV(0);
    scalar Tfavav(0);
    scalar Tfmaxav(0);
    scalar Tmavav(0);
    scalar TpSav(0);
    scalar TfSav(0);
    scalar Tmoutav(0);

    forAll(this->cellList_, i)
    {
        label celli(this->cellList_[i]);

        label regioni(cellToRegion_[celli]);
        //scalar fractionOfPowerFromNeutronics(fractionOfPowerFromNeutronics_[regioni]);
        scalar rcore(rcore_[regioni]);
        scalar rmatrix(rmatrix_[regioni]);
        scalar rshell(rshell_[regioni]);
        scalar rfuel(rfuel_[regioni]);
        scalar rbuffer(rbuffer_[regioni]);
        scalar rPyCin(rPyCin_[regioni]);
        scalar rSiC(rSiC_[regioni]);
        scalar rPyCout(rPyCout_[regioni]);

        label  Ntriso(Ntriso_[regioni]);

        scalar kfuel(kfuelCoeffs_[regioni].value(Tfav_[celli]));
        scalar kbuffer(kbufferCoeffs_[regioni].value((TfS_[celli]+Tmav_[celli])/2));
        scalar kSiC(kSiCCoeffs_[regioni].value((TfS_[celli]+Tmav_[celli])/2));
        scalar kPyC(kPyCCoeffs_[regioni].value((TfS_[celli]+Tmav_[celli])/2));
        scalar kgraphite(kgraphiteCoeffs_[regioni].value(Tmav_[celli]));

        scalar q = structure_.powerDensityNeutronics()[celli];
        const scalar& dV(V[celli]);
        scalar Q_cell = q*dV;
        scalar Q_pebble_tot =
            alpha_[celli]== 0 ?
            0
            :Q_cell*((4.0/3.0 * pi_ * pow(rshell, 3))/dV)*(1.0/alpha_[celli]);
        //scalar q_pebble = Q_pebble_tot / (4.0/3.0 * pi_ * pow(rshell, 3));
        scalar q_matrix = Q_pebble_tot / (4.0/3.0 * pi_ * (pow(rmatrix, 3) - pow(rcore, 3)));
        scalar Q_triso = Q_pebble_tot / Ntriso;
        scalar q_fuel = Q_triso / (4.0/3.0 * pi_ * pow(rfuel, 3));
        surfaceHeatFlux_[celli] = Q_pebble_tot / (4.0 * pi_ * pow(rshell, 2));

        // Calculation of effective conductivity in matrix (Maxwell model)
        scalar ktriso = (1/rfuel-1/rPyCout)/((1/kbuffer)*(1/rfuel-1/rbuffer)+(1/kPyC)*(1/rbuffer-1/rPyCin)+(1/kSiC)*(1/rPyCin-1/rSiC)+1/kPyC*(1/rSiC-1/rPyCout));
        scalar PF = (4.0/3.0*pi_*pow(rPyCout,3)) * static_cast<double >(Ntriso) / (4.0/3.0*pi_*(pow(rmatrix,3)-pow(rcore,3)));
        scalar kappa_coef = ktriso/kgraphite;
        scalar beta_coef  = (kappa_coef-1)/(kappa_coef+2);
        scalar keffmatrix = kgraphite*((1+2*beta_coef*PF)/(1-beta_coef*PF));

        // Calculation of average matrix and fuel temperatures directly
        TpS_[celli]   = (surfaceHeatFlux_[celli] + HTSum[celli]) / HSum[celli];
        Tmout_[celli] = TpS_[celli]   + (Q_pebble_tot/(4*pi_*kgraphite))*(1/rmatrix-1/rshell);
        Tmav_[celli]  = Tmout_[celli] + (q_matrix/(6*keffmatrix))*(pow(rmatrix,2)-3*(pow(rmatrix,5)-pow(rcore,5))/(5*(pow(rmatrix,3)-pow(rcore,3))))+(pow(rcore,3)*q_matrix)/(3*keffmatrix)*(1/rmatrix-3*(pow(rmatrix,2)-pow(rcore,2))/(2*(pow(rmatrix,3)-pow(rcore,3))));
        TfS_[celli]   = Tmav_[celli]  + Q_triso/(4*pi_) * (1/kbuffer*(1/rfuel-1/rbuffer) +  1/kPyC*(1/rbuffer-1/rPyCin) + 1/kSiC*(1/rPyCin-1/rSiC) + 1/kPyC*(1/rSiC-1/rPyCout));
        Tfav_[celli]  = TfS_[celli]   + q_fuel*pow(rfuel,2) / (15*kfuel);
        Tfmax_[celli] = TfS_[celli]   + q_fuel*pow(rfuel,2) / (6*kfuel);
        keffMatrix_[celli] = keffmatrix;

        rhofuel_[celli] = rhofuelCoeffs_[regioni].value(Tfav_[celli]);
        rhobuffer_[celli] = rhobufferCoeffs_[regioni].value((TfS_[celli]+Tmav_[celli])/2);
        rhoSiC_[celli] = rhoSiCCoeffs_[regioni].value((TfS_[celli]+Tmav_[celli])/2);
        rhoPyC_[celli] = rhoPyCCoeffs_[regioni].value((TfS_[celli]+Tmav_[celli])/2);
        rhographite_[celli] = rhographiteCoeffs_[regioni].value(Tmav_[celli]);
        rhocoolant_[celli] = rhocoolantCoeffs_[regioni].value(HTSum[celli] / HSum[celli]);

        totV += dV;
        Tfavav += Tfav_[celli]*dV;
        Tfmaxav += Tfmax_[celli]*dV;
        Tmavav += Tmav_[celli]*dV;
        TpSav += TpS_[celli]*dV;
        TfSav += TfS_[celli]*dV;
        Tmoutav += Tmout_[celli]*dV;

        // Update overall minimums and maximums
        updateLocalAvgGlobalMinMaxT(Tfav_[celli],  Tfavmin_,  Tfavmax_);
        updateLocalAvgGlobalMinMaxT(Tfmax_[celli],  Tfmaxmin_,  Tfmaxmax_);
        updateLocalAvgGlobalMinMaxT(Tmout_[celli], Tmoutmin_, Tmoutmax_);
        updateLocalAvgGlobalMinMaxT(Tmav_[celli],  Tmavmin_,  Tmavmax_);
        updateLocalAvgGlobalMinMaxT(TpS_[celli],   TpSmin_,   TpSmax_);
        updateLocalAvgGlobalMinMaxT(TfS_[celli],   TfSmin_,   TfSmax_);
    }

    //- Save these to the dictionary
    this->IOdictionary::set("Tfav", Tfav_);
    this->IOdictionary::set("Tfmax", Tfmax_);
    this->IOdictionary::set("Tmav", Tmav_);
    this->IOdictionary::set("TpS", TpS_);
    this->IOdictionary::set("Tmout", Tmout_);
    this->IOdictionary::set("TfS", TfS_);
    this->IOdictionary::set("keffmatrix", keffMatrix_);
    this->IOdictionary::set("surfaceHeatFlux", surfaceHeatFlux_);
    this->IOdictionary::set("rhofuel", rhofuel_);
    this->IOdictionary::set("rhobuffer", rhobuffer_);
    this->IOdictionary::set("rhoSiC", rhoSiC_);
    this->IOdictionary::set("rhoPyC", rhoPyC_);
    this->IOdictionary::set("rhographite", rhographite_);
    this->IOdictionary::set("rhocoolant", rhocoolant_);

    reduce(totV, sumOp<scalar>());
    reduce(Tfavav, sumOp<scalar>());
    reduce(Tfmaxav, sumOp<scalar>());
    reduce(Tmavav, sumOp<scalar>());
    reduce(TpSav, sumOp<scalar>());
    reduce(Tmoutav, sumOp<scalar>());
    reduce(TfSav, sumOp<scalar>());
    Tfavav /= totV;
    Tfmaxav /= totV;
    Tmavav /= totV;
    TpSav /= totV;
    Tmoutav /= totV;
    TfSav /= totV;

    reduce(Tfavmax_,  maxOp<scalar>());
    reduce(Tfavmin_,  minOp<scalar>());
    reduce(Tfmaxmax_,  maxOp<scalar>());
    reduce(Tfmaxmin_,  minOp<scalar>());
    reduce(Tmoutmax_, maxOp<scalar>());
    reduce(Tmoutmin_, minOp<scalar>());
    reduce(Tmavmax_,  maxOp<scalar>());
    reduce(Tmavmin_,  minOp<scalar>());
    reduce(TpSmax_,   maxOp<scalar>());
    reduce(TpSmin_,   minOp<scalar>());
    reduce(TfSmax_,   maxOp<scalar>());
    reduce(TfSmin_,   minOp<scalar>());

    Info<< "T.nuclearSteadyStatePebble.TpS (avg min max) = "    << TpSav    << " " << TpSmin_    << " " << TpSmax_    << " K" << endl;
    Info<< "T.nuclearSteadyStatePebble.Tmout (avg min max) = "  << Tmoutav  << " " << Tmoutmin_  << " " << Tmoutmax_  << " K" << endl;
    Info<< "T.nuclearSteadyStatePebble.Tmav (avg min max) = "   << Tmavav   << " " << Tmavmin_   << " " << Tmavmax_   << " K" << endl;
    Info<< "T.nuclearSteadyStatePebble.TfS (avg min max) = "    << TfSav    << " " << TfSmin_    << " " << TfSmax_    << " K" << endl;
    Info<< "T.nuclearSteadyStatePebble.Tfav (avg min max) = "   << Tfavav   << " " << Tfavmin_   << " " << Tfavmax_   << " K" << endl;
    Info<< "T.nuclearSteadyStatePebble.Tfmax (avg min max) = "  << Tfmaxav  << " " << Tfmaxmin_  << " " << Tfmaxmax_  << " K" << endl;

    //- Save these to the dictionary
    this->IOdictionary::set("Tfavav",     Tfavav);
    this->IOdictionary::set("Tfmaxav",    Tfmaxav);
    this->IOdictionary::set("Tmavav",     Tmavav);
    this->IOdictionary::set("TpSav",      TpSav);
    this->IOdictionary::set("Tmoutav",    Tmoutav);
    this->IOdictionary::set("TfSav",      TfSav);
    this->IOdictionary::set("Tfavmax",    Tfavmax_);
    this->IOdictionary::set("Tfavmin",    Tfavmin_);
    this->IOdictionary::set("Tfmaxmax",   Tfavmax_);
    this->IOdictionary::set("Tfmaxmin",   Tfavmin_);
    this->IOdictionary::set("Tmoutmax",   Tmoutmax_);
    this->IOdictionary::set("Tmoutmin",   Tmoutmin_);
    this->IOdictionary::set("Tmavmax",    Tmavmax_);
    this->IOdictionary::set("Tmavmin",    Tmavmin_);
    this->IOdictionary::set("TpSmax",     TpSmax_);
    this->IOdictionary::set("TpSmin",     TpSmin_);
    this->IOdictionary::set("TfSmax",     TfSmax_);
    this->IOdictionary::set("TfSmin",     TfSmin_);
}

void Foam::powerModels::nuclearSteadyStatePebble::correctT(volScalarField& T) const
{
    //- Set T to pebble surface temperature, i.e. TpS_
    forAll(cellList_, i)
    {
        label celli(cellList_[i]);
        T[celli] = TpS_[celli];
    }
}

void Foam::powerModels::nuclearSteadyStatePebble::updateLocalAvgGlobalMinMaxT
(
    scalar& T,
    scalar& Tmin,
    scalar& Tmax
)
{
    if (T > Tmax) Tmax = T;
    if (T < Tmin) Tmin = T;
}
// ************************************************************************* //
