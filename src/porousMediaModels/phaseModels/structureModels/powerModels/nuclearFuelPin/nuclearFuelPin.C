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

#include "nuclearFuelPin.H"
#include "structure.H"
#include "addToRunTimeSelectionTable.H"
#include "SquareMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace powerModels
{
    defineTypeNameAndDebug(nuclearFuelPin, 0);
    addToRunTimeSelectionTable
    (
        powerModel,
        nuclearFuelPin,
        powerModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::powerModels::nuclearFuelPin::nuclearFuelPin
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
    Trad_
    (
        IOobject
        (
            "Trad."+typeName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_.cells().size()
    ),
    Tfi_
    (
        IOobject
        (
            "Tfi."+typeName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    Tfo_
    (
        IOobject
        (
            "Tfo."+typeName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    Tci_
    (
        IOobject
        (
            "Tci."+typeName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    Tco_
    (
        IOobject
        (
            "Tco."+typeName,
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
    Tcav_
    (
        IOobject
        (
            "Tcav."+typeName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    Tfmax_(0),
    Tfmin_(1e69),
    Tcmax_(0),
    Tcmin_(1e69),
    fractionOfPowerFromNeutronics_(0),
    fuelMeshSize_(0),
    cladMeshSize_(0),
    meshSize_(0),
    r_(0),
    rfi_(0),
    rfo_(0),
    rci_(0),
    rco_(0),
    drf_(0),
    drc_(0),
    drg_(0),
    rhoCpf_(0),
    rhoCpc_(0),
    kf_(0),
    kc_(0),
    gapH_(0),
    hollowFuel_(0),
    cellToRegion_(mesh_.cells().size(), 0),
    regionIndexToRegionName_(0),
    pi_(constant::mathematical::pi),
    dA_(0),
    gapHPowerDensityTable_(0),
    useGapHPowerDensityTable_(0)
{
    structure_.setRegionField(*this, structureRef.powerDensityNeutronics(), "powerDensity");

    bool foundBoundaryTemperatures
    (
           Tfi_.typeHeaderOk<volScalarField>(true)
        && Tfo_.typeHeaderOk<volScalarField>(true)
        && Tci_.typeHeaderOk<volScalarField>(true)
        && Tco_.typeHeaderOk<volScalarField>(true)
    );

    typedef IOFieldField<Field, scalar> scalarFieldField;
    bool foundTrad(Trad_.typeHeaderOk<scalarFieldField>(true));

    scalarList Tf0(0);
    scalarList Tc0(0);

    forAll(this->toc(), regioni)
    {
        word region(this->toc()[regioni]);
        const dictionary& dict(this->subDict(region));

        // Setup cellToRegion_ mapping
        const labelList& regionCells
        (
            structure_.cellLists()[region]
        );
        forAll(regionCells, i)
        {
            label celli(regionCells[i]);
            cellToRegion_[celli] = regioni;
        }

        // Add to regionIndexToRegionName_ mapping
        regionIndexToRegionName_.append(region);

        // Read region dict entries
        scalar fractionOfPowerFromNeutronics(dict.lookupOrDefault<scalar>("fractionOfPowerFromNeutronics",1.0));
        scalar rfi(dict.get<scalar>("fuelInnerRadius"));
        scalar rfo(dict.get<scalar>("fuelOuterRadius"));
        scalar rci(dict.get<scalar>("cladInnerRadius"));
        scalar rco(dict.get<scalar>("cladOuterRadius"));
        scalar fuelFraction((pow(rfo,2)-pow(rfi,2)) / (pow(rco,2)));
        label fuelMeshSize(dict.get<label>("fuelMeshSize"));
        label cladMeshSize(dict.get<label>("cladMeshSize"));
        label meshSize(fuelMeshSize + cladMeshSize);
        scalar drf((rfo-rfi) / (fuelMeshSize-1));
        scalar drc((rco-rci) / (cladMeshSize-1));
        scalar drg(rci - rfo);
        scalar rhoCpf(0);
        if (dict.found("fuelRho") && dict.found("fuelCp"))
        {
            rhoCpf =
                dict.get<scalar>("fuelRho")*
                dict.get<scalar>("fuelCp");
        }
        else if (dict.found("fuelRhoCp"))
        {
            rhoCpf = dict.get<scalar>("fuelRhoCp");
        }
        else
        {
            FatalErrorInFunction
                << "nuclearFuelPin region: " << region << " -> "
                << "specify either fuelRhoCp or both fuelRho and fuelCp"
                << exit(FatalError);
        }

        scalar rhoCpc(0);
        if (dict.found("cladRho") and dict.found("cladCp"))
        {
            rhoCpc =
                dict.get<scalar>("cladRho")*
                dict.get<scalar>("cladCp");
        }
        else if (dict.found("cladRhoCp"))
        {
            rhoCpc = dict.get<scalar>("cladRhoCp");
        }
        else
        {
            FatalErrorInFunction
                << "nuclearFuelPin region: " << region << " -> "
                << "specify either cladRhoCp or both cladRho and cladCp"
                << exit(FatalError);
        }
        scalar kf(dict.get<scalar>("fuelK"));
        scalar kc(dict.get<scalar>("cladK"));

        bool hollowFuel((rfi >= 1e-5) ? true : false);

        if (!foundBoundaryTemperatures && !foundTrad)
        {
            Tf0.append(dict.get<scalar>("fuelT"));
            Tc0.append(dict.get<scalar>("cladT"));
        }

        // Calc mesh array
        scalarList r(0);
        r.append(rfi);
        for (int i = 0; i < fuelMeshSize-1; i++)
        {
            r.append(r.last() + drf);
        }
        r.append(r.last() + drg);
        for (int i = 0; i < cladMeshSize-1; i++)
        {
            r.append(r.last() + drc);
        }

        // Calc dA
        scalarList dA(0);
        dA.append(pi_*(sqr(r[0]+drf/2.0)-sqr(r[0])));
        for(int i = 1; i < fuelMeshSize-1; i++)
        {
            dA.append(pi_*(sqr(r[i]+drf/2.0)-sqr(r[i]-drf/2.0)));
        }
        dA.append
        (
            pi_*(sqr(r[fuelMeshSize-1])-sqr(r[fuelMeshSize-1]-drf/2.0))
        );
        dA.append
        (
            pi_*(sqr(r[fuelMeshSize]+drf/2.0)-sqr(r[fuelMeshSize]))
        );
        for(int i = 1; i < cladMeshSize-1; i++)
        {
            dA.append(pi_*(sqr(r[i]+drc/2.0)-sqr(r[i]-drc/2.0)));
        }
        dA.append
        (
            pi_*(sqr(r[meshSize-1])-sqr(r[meshSize-1]-drc/2.0))
        );

        // Fill in lists for this region
        fractionOfPowerFromNeutronics_.append(fractionOfPowerFromNeutronics),
        fuelMeshSize_.append(fuelMeshSize);
        cladMeshSize_.append(cladMeshSize);
        meshSize_.append(meshSize);
        r_.append(r);
        rfi_.append(rfi);
        rfo_.append(rfo);
        rci_.append(rci);
        rco_.append(rco);

        drf_.append(drf);
        drc_.append(drc);
        drg_.append(drg);
        rhoCpf_.append(rhoCpf);
        rhoCpc_.append(rhoCpc);
        kf_.append(kf);
        kc_.append(kc);

        hollowFuel_.append(hollowFuel);
        dA_.append(dA);

        // Construct gapHPowerDensityTable if found, otherwise use the
        // provided constant value
        scalar gapH(0);
        word tableName("gapHPowerDensity");
        bool foundTable(dict.found(tableName));
        bool foundValue(dict.found("gapH"));
        if (foundTable)
        {
            /*
                The table Function1 requires a particular input. In general,
                it should be something like this

                topLevelDictName
                {
                    type                table;
                    topLevelDictName    table
                    (
                        (0 0)
                        (1 1)
                        (...)
                    );
                }

                however, I only want the user to give an input in the form

                powerModel
                {
                    type        nuclearFuelPin;

                    ...

                    gapHPowerDensity table
                    (
                        (0 0)
                        (1 1)
                        (...)
                    );
                }

                in which powerModel is the topLevelDictName.

                The code below does just this, by creating a copy dict of
                powerModel renaming it to gapHPowerDensity table, resetting
                type to table, and passing that to the Function1 table
                selector
            */
            dictionary tableDict(tableName);
            tableDict.merge(dict);
            tableDict.set("type", "table");
            gapHPowerDensityTable_.insert
            (
                region,
                autoPtr<Function1<scalar>>
                (
                    Function1<scalar>::New
                    (
                        tableName,
                        tableDict,
                        "table"
                    )
                )
            );
        }
        if (foundValue)
        {
            gapH = dict.get<scalar>("gapH");
        }

        // The gapH list needs to have the same length as the number of
        // regions no matter what, or the indexing will stop working as
        // intended
        gapH_.append(gapH);

        useGapHPowerDensityTable_.append(foundTable);

        fuelFraction_.append(fuelFraction);

        if (foundTable && foundValue)
        {
            FatalErrorInFunction
                << "nuclearFuelPin region: " << region << " -> "
                << "provide either a gapH value or a gapHPowerDensityTable "
                << "but not both!"
                << exit(FatalError);
        }
    }

    // If Trad not found, init it from either the boundary temperatures
    // (I mean boundary in a mathematical sense, i.e. inner/outer fuel/clad
    // temperature) or from dictionary values Tf0, Tc0 read previously
    if (!foundTrad)
    {
        forAll(mesh_.cells(), i)
        {
            Trad_.set(i, new Field<scalar>(0, 0));
        }

        // If the files are present, reconstruct initial Trad_ profile
        // analytically. The analytical form is:
        // T(r) = -(1/4)*powerDensity_(r)*r^2/k + C*ln(r)/k + D
        // with C and D coming from imposing fixedValue BC on all sides,
        // equal to the starting temperatures found in the files.
        if (foundBoundaryTemperatures)
        {
            Info<< "Found nuclearFuelPin temperatures reconstructing profiles"
                << endl;
            forAll(this->cellList_, i)
            {
                label celli(this->cellList_[i]);

                label regioni(cellToRegion_[celli]);
                const scalarList& rList(r_[regioni]);
                scalar kf(kf_[regioni]);
                scalar kc(kc_[regioni]);
                scalar rfi(rfi_[regioni]);
                scalar rfo(rfo_[regioni]);
                scalar rci(rci_[regioni]);
                scalar rco(rco_[regioni]);
                scalar fractionOfPowerFromNeutronics(fractionOfPowerFromNeutronics_[regioni]);
                bool hollowFuel(hollowFuel_[regioni]);

                Trad_.set(celli, new Field<scalar>(meshSize_[regioni], 0));
                scalar q = structure_.powerDensityNeutronics()[celli] * fractionOfPowerFromNeutronics;
                scalar tfi = Tfi_[celli];
                scalar tfo = Tfo_[celli];
                scalar tci = Tci_[celli];
                scalar tco = Tco_[celli];
                scalar Cf;
                scalar Df;
                scalar Cc;
                scalar Dc;

                if (hollowFuel)
                {
                    Cf =
                        (kf*(tfi-tfo)-0.25*q*(sqr(rfo)-sqr(rfi)))/
                        log(rfi/rfo);
                }
                else
                {
                    Cf = 0.0;
                }
                Df = tfo + (0.25*q*sqr(rfo)-Cf*log(rfo))/kf;
                Cc = (tco-tci)*kc/(log(rco/rci));
                Dc = tci - log(rci)*Cc/kc;

                forAll(Trad_[celli], j)
                {
                    scalar r(rList[j]);
                    if (j < fuelMeshSize_[regioni])
                    {
                        Trad_[celli][j] = -0.25*q*sqr(r)/kf + Df;
                        Trad_[celli][j] += (hollowFuel)
                            ? Cf*log(r)/kf
                            : 0.0;
                    }
                    else
                    {
                        Trad_[celli][j] = Cc*log(r)/kc + Dc;
                    }
                }
           }
        }
        else // Otherwise, read from dict
        {
            Info<< "Reading nuclearFuelPin initial temperatures from dictionary"
                << endl;

            forAll(this->cellList_, i)
            {
                label celli(this->cellList_[i]);
                label regioni(cellToRegion_[celli]);
                Trad_.set(celli, new Field<scalar>(meshSize_[regioni], 0));
                forAll(Trad_[celli], subCelli)
                {
                    Trad_[celli][subCelli] = (subCelli < fuelMeshSize_[regioni])
                        ? Tf0[regioni]
                        : Tc0[regioni];
                }
            }
        }
    }
    else
    {
        Info<< "Setting nuclearFuelPin initial temperatures from "
            << Trad_.name()
            << endl;
    }

    // Set I/O fields and compute initial scalar max, min
    scalar Tfavav(0);
    scalar Tcavav(1e69);
    scalar totV(0.0);
    const scalarList& V(mesh_.V());
    forAll(this->cellList_, i)
    {
        label celli(this->cellList_[i]);
        label regioni(cellToRegion_[celli]);
        const label& nf(fuelMeshSize_[regioni]);
        const label& n(meshSize_[regioni]);
        const scalarField& Trad(Trad_[celli]);
        Tfi_[celli] = Trad[0];
        Tfo_[celli] = Trad[nf-1];
        Tci_[celli] = Trad[nf];
        Tco_[celli] = Trad[n-1];

        const scalarList& rRegion(r_[regioni]);
        scalar& Tfavi(Tfav_[celli]);
        scalar& Tcavi(Tcav_[celli]);

        updateLocalAvgGlobalMinMaxT
        (
            0,
            nf,
            rRegion,
            drf_[regioni],
            Trad,
            Tfavi,
            Tfmin_,
            Tfmax_
        );
        updateLocalAvgGlobalMinMaxT
        (
            nf,
            n,
            rRegion,
            drc_[regioni],
            Trad,
            Tcavi,
            Tcmin_,
            Tcmax_
        );

        // Update average fuel and clad temp used for coupling
        this->structureRef().TFuelAv()[celli] = Tfav_[celli];
        this->structureRef().TCladAv()[celli] = Tcav_[celli];

        // This is for updating the global averages, not the local cell ones!
        const scalar& dV(V[celli]);
        totV += dV;
        Tfavav += Tfavi*dV;
        Tcavav += Tcavi*dV;
    }

    // Sync across processors
    reduce(totV, sumOp<scalar>());
    reduce(Tfavav, sumOp<scalar>());
    reduce(Tcavav, sumOp<scalar>());
    reduce(Tfmax_, maxOp<scalar>());
    reduce(Tfmin_, minOp<scalar>());
    reduce(Tcmax_, maxOp<scalar>());
    reduce(Tcmin_, minOp<scalar>());

    Tfavav /= totV;
    Tcavav /= totV;

    // Initialize in dict
    this->IOdictionary::set("Tfavav", Tfavav);
    this->IOdictionary::set("Tcavav", Tcavav);
    this->IOdictionary::set("Tfmax", Tfmax_);
    this->IOdictionary::set("Tfmin", Tfmin_);
    this->IOdictionary::set("Tcmax", Tcmax_);
    this->IOdictionary::set("Tcmin", Tcmin_);

    Tfi_.correctBoundaryConditions();
    Tfo_.correctBoundaryConditions();
    Tci_.correctBoundaryConditions();
    Tco_.correctBoundaryConditions();

    // Finally, set up interfacial area
    this->setInterfacialArea();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::powerModels::nuclearFuelPin::~nuclearFuelPin()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::powerModels::nuclearFuelPin::setInterfacialArea()
{
    forAll(this->cellList_, i)
    {
        const label& celli(this->cellList_[i]);
        iA_[celli] = 2.0*alpha_[celli]/(rco_[cellToRegion_[celli]]);
    }
    iA_.correctBoundaryConditions();
}

void Foam::powerModels::nuclearFuelPin::updateLocalAvgGlobalMinMaxT
(
    const label& starti,
    const label& endi,
    const scalarList& r,
    const scalar& dr,
    const scalarField& Trad,
    scalar& Tavi,
    scalar& Tmin,
    scalar& Tmax
)
{
    scalar intr(0);
    scalar intTr(0);
    for (int j = starti; j < endi; j++)
    {
        const scalar& T(Trad[j]);
        scalar rdr(r[j]*dr);

        // Cells at the mesh ends are only half as wide (the other half
        // belongs to the ghost node). Thus, weigh temperatures at the extrema
        // by a factor 0.5
        if (j == starti or j == endi-1)
        {
            intr += rdr / 2.0;
            intTr += T*rdr / 2.0;
        }
        else
        {
            intr += rdr;
            intTr += T*rdr;
        }
        if (T > Tmax) Tmax = T;
        if (T < Tmin) Tmin = T;
    }
    Tavi = intTr/intr;
}

void
Foam::powerModels::nuclearFuelPin::updateLocalTemperatureProfile
(
    const label& celli,
    const scalar& HTSumi,
    const scalar& HSumi
)
{
    scalarField& Trad(Trad_[celli]);

    // Read region values
    const label& regioni(cellToRegion_[celli]);
    const word& region(regionIndexToRegionName_[regioni]);
    const label& fuelMeshSize(fuelMeshSize_[regioni]);
    const label& meshSize(meshSize_[regioni]);
    const scalarList& rRegion(r_[regioni]);
    const scalarList& dARegion(dA_[regioni]);
    const scalar& drf(drf_[regioni]);
    const scalar& drc(drc_[regioni]);
    const scalar& kf(kf_[regioni]);
    const scalar& kc(kc_[regioni]);
    const scalar& rfo(rfo_[regioni]);
    const scalar& rci(rci_[regioni]);
    const scalar& fractionOfPowerFromNeutronics(fractionOfPowerFromNeutronics_[regioni]);

    const scalarField& TOld = Trad_.oldTime()[celli];

    // Update power density
    const scalar& qRef(structure_.powerDensityNeutronics()[celli]);
    scalar q = qRef * fractionOfPowerFromNeutronics;

    scalar gapH
    (
        (useGapHPowerDensityTable_[regioni])
            ? gapHPowerDensityTable_[region]->value(q)
            : gapH_[regioni]
    );

    // Init matrix, source
    SquareMatrix<scalar> M(meshSize, meshSize, Foam::zero());
    List<scalar> S(meshSize, 0.0);

    // Recurrent quantities
    scalar dt(mesh_.time().deltaT().value());
    scalar Xf(rhoCpf_[regioni]/dt);
    scalar Xc(rhoCpc_[regioni]/dt);
    scalar twoPkByDrf(2.0*pi_*kf/drf);
    scalar twoPkByDrc(2.0*pi_*kc/drc);
    scalar drhf(drf/2.0);
    scalar drhc(drc/2.0);

    // Construct matrix, source
    {
        // Set zeroGradient BC at fuel inner surface
        {
            const scalar& r(rRegion[0]);
            const scalar& dA(dARegion[0]);
            scalar B(twoPkByDrf*(r+drhf));
            scalar XdA(Xf*dA);
            M[0][1] = -B;
            M[0][0] = B+XdA;
            S[0] =    q*dA+TOld[0]*XdA;
        }

        // Fuel bulk
        for (int i = 1; i < fuelMeshSize-1; i++)
        {
            const scalar& r(rRegion[i]);
            const scalar& dA(dARegion[i]);
            scalar B(twoPkByDrf*(r+drhf));
            scalar C(twoPkByDrf*(r-drhf));
            scalar XdA(Xf*dA);
            M[i][i+1] = -B;
            M[i][i-1] = -C;
            M[i][i] =   B+C+XdA;
            S[i] =      q*dA+TOld[i]*XdA;
        }

        // Fuel outer surface, convective BC with inner cladding surface via
        // gap conductance
        {
            label i(fuelMeshSize-1);
            const scalar& r(rRegion[i]);
            const scalar& dA(dARegion[i]);
            scalar C(twoPkByDrf*(r-drhf));
            scalar D(2.0*pi_*r*gapH);
            scalar XdA(Xf*dA);
            M[i][i+1] = -D;
            M[i][i-1] = -C;
            M[i][i] =   C+D+XdA;
            S[i] =      q*dA+TOld[i]*XdA;
        }

        // Cladding inner surface, convective BC with outer fuel surface via
        // gap conductance adjusted by radii ratio (to preserve total heat
        // flow as geometry is cylindrical)
        {
            label i(fuelMeshSize);
            const scalar& r(rRegion[i]);
            const scalar& dA(dARegion[i]);
            scalar C(twoPkByDrc*(r+drhc));
            scalar D(2.0*pi_*r*(rfo/rci)*gapH);
            scalar XdA(Xc*dA);
            M[i][i+1] = -C;
            M[i][i-1] = -D;
            M[i][i] =   C+D+XdA;
            S[i] =      TOld[i]*XdA;
        }

        // Cladding bulk
        for (int i = fuelMeshSize+1; i < meshSize-1; i++)
        {
            const scalar& r(rRegion[i]);
            const scalar& dA(dARegion[i]);
            scalar B(twoPkByDrc*(r+drhc));
            scalar C(twoPkByDrc*(r-drhc));
            scalar XdA(Xc*dA);
            M[i][i+1] = -B;
            M[i][i-1] = -C;
            M[i][i] =   B+C+XdA;
            S[i] =      TOld[i]*XdA;
        }

        // Cladding outer surface, convective BC with fluid(s) wetting the pin
        {
            label i(meshSize-1);
            const scalar& r(rRegion[i]);
            const scalar& dA(dARegion[i]);
            scalar C(twoPkByDrc*(r-drhc));
            scalar D(2.0*pi_*r);
            scalar XdA(Xc*dA);
            M[i][i-1] = -C;
            M[i][i] =   C+D*HSumi+XdA;
            S[i] =      TOld[i]*XdA + D*HTSumi;
        }
    }

    // Solve linear system
    solve(Trad, M, S);

    // Set fields (inner/outer fuel/clad)
    Tfi_[celli] = Trad[0];
    Tfo_[celli] = Trad[fuelMeshSize-1];
    Tci_[celli] = Trad[fuelMeshSize];
    Tco_[celli] = Trad[meshSize-1];

    // Update local T averages and local min/max
    scalar& Tfavi(Tfav_[celli]);
    scalar& Tcavi(Tcav_[celli]);
    updateLocalAvgGlobalMinMaxT
    (
        0,
        fuelMeshSize,
        rRegion,
        drf,
        Trad,
        Tfavi,
        Tfmin_,
        Tfmax_
    );
    updateLocalAvgGlobalMinMaxT
    (
        fuelMeshSize,
        meshSize,
        rRegion,
        drc,
        Trad,
        Tcavi,
        Tcmin_,
        Tcmax_
    );

    // Update average fuel and clad temp used for coupling
    this->structureRef().TFuelAv()[celli] = Tfav_[celli];
    this->structureRef().TCladAv()[celli] = Tcav_[celli];

    /*
    // Check energy conservation via linear power comparison (analytic
    // vs numerical at cladding surface)
    const scalar& rfi(rfi_[regioni]);
    const scalar& rco(rco_[regioni]);
    scalar analyticalLP
    (
        q*pi_*(sqr(rfo)-sqr(rfi))
    );
    scalar numericalLP
    (
        (HSumi*Tco_[celli]-HTSumi)*2.0*pi_*rco
    );
    Info<< celli << " " << numericalLP << " " << analyticalLP << " W/m"
        << endl;

    // Check energy conservation via heat flux comparison. However, since
    // you can't directly compare heat fluxes (it's total heat that is
    // conserved, not heat fluxes), the heat flux trhough the cladding is
    // adjusted to take into consideration the difference in surface areas
    // between outer cladding and outer fuel
    scalar heatFluxg(gapH*(Tfo_[celli]-Tci_[celli]));
    scalar heatFluxAdjc((HSumi*Tco_[celli]-HTSumi)*(rco/rfo));
    Info<< celli << " " << heatFluxg << " " << heatFluxAdjc << " W/m2"
        << endl;
    */
}


void Foam::powerModels::nuclearFuelPin::correct
(
    const volScalarField& HTSum,  // == SUM_j [htc_j*T_j*frac_j]
    const volScalarField& HSum    // == SUM_j [htc_j*frac_j]
)
{

    forAll(this->toc(), regioni)
    {
        if(powerOffCriterionModelPtr_.set(regioni))
        {
            if(powerOffCriterionModelPtr_[regioni].powerOffCriterion())
            {
                word region(this->toc()[regioni]);
                const labelList& regionCells
                (
                    structure_.cellLists()[region]
                );
                forAll(regionCells, j)
                {
                    label cellj(regionCells[j]);
                    structure_.powerDensityNeutronics()[cellj] = 0.0;
                }
            }
        }
    }
    // Reset min, max, fuel, clad temperatures
    Tfmax_ = 0.0;
    Tfmin_ = 1e69;
    Tcmax_ = 0.0;
    Tcmin_ = 1e69;
    scalar totalPower(0.0);

    // Update temperatures cell-by-cell and compute averages over the entire
    // spatial extent of the nuclearFuelPin model (what I call global
    // averages, opposed to local averages, which are the average temperature
    // values, fuel and clad, of the local radial pin temperature profile)
    const scalarField& V(mesh_.V());
    scalar totV(0);
    scalar Tfavav(0);
    scalar Tcavav(0);
    forAll(this->cellList_, i)
    {
        label celli(this->cellList_[i]);
        label regioni(cellToRegion_[celli]);
        updateLocalTemperatureProfile(celli, HTSum[celli], HSum[celli]);
        const scalar& dV(V[celli]);
        totV += dV;
        Tfavav += Tfav_[celli]*dV;
        Tcavav += Tcav_[celli]*dV;
        totalPower += structure_.powerDensityNeutronics()[celli] * alpha_[celli] * V[celli]*fuelFraction_[regioni];
    }
    reduce(totV, sumOp<scalar>());
    reduce(Tfavav, sumOp<scalar>());
    reduce(Tcavav, sumOp<scalar>());
    reduce(totalPower, sumOp<scalar>());
    Tfavav /= totV;
    Tcavav /= totV;

    reduce(Tfmax_, maxOp<scalar>());
    reduce(Tfmin_, minOp<scalar>());
    reduce(Tcmax_, maxOp<scalar>());
    reduce(Tcmin_, minOp<scalar>());

    Info<< "T.nuclearFuelPin.fuel (avg min max) = "
        << Tfavav << " " << Tfmin_ << " " << Tfmax_ << " K"
        << endl;
    Info<< "T.nuclearFuelPin.clad (avg min max) = "
        << Tcavav << " " << Tcmin_ << " " << Tcmax_ << " K"
        << endl;
    Info<< "Total power in nuclearFuelPin = " << totalPower << " W"
        << endl;

    // Save these to the dictionary
    this->IOdictionary::set("Tfavav", Tfavav);
    this->IOdictionary::set("Tcavav", Tcavav);
    this->IOdictionary::set("Tfmax", Tfmax_);
    this->IOdictionary::set("Tfmin", Tfmin_);
    this->IOdictionary::set("Tcmax", Tcmax_);
    this->IOdictionary::set("Tcmin", Tcmin_);
}


void Foam::powerModels::nuclearFuelPin::correctT(volScalarField& T) const
{
    // Set T to pin surface temperature, i.e. Tco_
    forAll(cellList_, i)
    {
        label celli(cellList_[i]);
        T[celli] = Tco_[celli];
    }
}


// ************************************************************************* //
