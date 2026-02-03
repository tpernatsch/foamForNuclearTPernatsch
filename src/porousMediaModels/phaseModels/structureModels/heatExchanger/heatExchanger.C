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

#include "heatExchanger.H"

// From forward declarations
#include "fluid.H"
#include "structure.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(heatExchanger, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatExchanger::heatExchanger
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    IOdictionary
    (
        IOobject
        (
            typeName,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dict
    ),
    mesh_(mesh),
    iA_(this->get<scalar>("volumetricArea")),
    Hw_(max(this->get<scalar>("wallConductance"), 1e-69)),
    primaryCells_( mesh_.cellZones()[this->get<word>("primary")]),
    secondaryCells_(0),
    constructed_(false),
    thisDictionary_(dict)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::heatExchanger::iA() const
{
    tmp<volScalarField> tiA
    (
        new volScalarField
        (
            IOobject
            (
                "",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("", dimArea/dimVol, 0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volScalarField& iA(tiA.ref());

    forAll(primaryCells_, i)
    {
        const label& celli(primaryCells_[i]);
        iA[celli] = iA_;
    }

    iA.correctBoundaryConditions();

    return tiA;
}


void Foam::heatExchanger::correct
(
    const volScalarField& HT,
    const volScalarField& H,
    volScalarField& THX
)
{
    // Get secondary side heat transfer information
    if (!constructed_)
    {
        constructSecondaryHX(constructed_);
    }

    const fvMesh& sMesh_ = mesh_.time().lookupObjectRef<fvMesh>(this->getOrDefault<word>("secondaryRegion", mesh_.name()));

    volScalarField sH
    (
        sMesh_.foundObject<volScalarField>("htc")
            ? sMesh_.lookupObject<volScalarField>("htc")
            : sMesh_.lookupObject<volScalarField>("htc.liquid.structure")
    );
    volScalarField sHT
    (
        sMesh_.foundObject<volScalarField>("htc")
            ? sMesh_.lookupObject<volScalarField>("htc")*sMesh_.lookupObject<volScalarField>("T")
            : sMesh_.lookupObject<volScalarField>("htc.liquid.structure")*sMesh_.lookupObject<volScalarField>("T.liquid")
    );

    // If twoPhase, correct the herat transfer parameters to account for both
    // phases

    if (sMesh_.foundObject<volScalarField>("htc.vapour.structure"))
    {
        const volScalarField& sT2 = (sMesh_.lookupObject<volScalarField>("T.vapour"));
        const volScalarField& sH2=(sMesh_.lookupObject<volScalarField>("htc.vapour.structure"));
        sHT += sH2 * sT2;
        sH += sH2;
    }

    scalarField HTpOs=mappingPtr_->mapTgtToSrc(HT.internalField());
    scalarField HTsOp = mappingPtr_->mapSrcToTgt(sHT.internalField());
    scalarField HpOs = mappingPtr_->mapTgtToSrc(H.internalField());
    scalarField HsOp = mappingPtr_->mapSrcToTgt(sH.internalField());

    // The max against 1e-69 is just to avoid the case in which the
    // heat transfer coefficient on both side is 0 (which results in
    // Ap*As = 1.0), a case in which no heat is to be transferred via the HX
    // but that would result in a division by 0 if not accounted for. Please
    // note that in such a scenario THX is set to 0 yet it does not matter
    // as said earlier, this only happens where the heat transfer coefficient
    // 0 or very close to it
    forAll(primaryCells_, i)
    {
        const label& celli(primaryCells_[i]);
        scalar Ap((Hw_+H[celli])/Hw_);
        scalar Bp(HT[celli]/Hw_);
        scalar As((Hw_+HsOp[celli])/Hw_);
        scalar Bs(HTsOp[celli]/Hw_);
        THX[celli] = (Bp*As+Bs)/max(Ap*As-1.0, 1e-69);
    }

    THX.correctBoundaryConditions();
}

void Foam::heatExchanger::constructSecondaryHX(bool& constructed)
{
    Info<< "Creating heatExchanger: " << thisDictionary_.dictName() << endl;

    const fvMesh& sMesh_ = mesh_.time().lookupObjectRef<dynamicFvMesh>(this->getOrDefault<word>("secondaryRegion", mesh_.name()));
    // Read secondary cellZones and their cells
    secondaryCells_ = sMesh_.cellZones()[this->get<word>("secondary")];

    // Calculate displacement vecdtor that, if applied to the primary
    // cellZone, would translate it to the secondary cellZone. This assumes
    // that the cellZones have the same shape and volume. If that is not
    // the case, a check is done afterwards (currently only on the volume, not
    // on the shape, which I will implement via comparing bounding boxes at
    // some point)
    vector pCOV(vector::zero);
    vector sCOV(vector::zero);
    const vectorField& primC(mesh_.C());
    const scalarField& primV(mesh_.V());
    const vectorField& secC(sMesh_.C());
    const scalarField& secV(sMesh_.V());
    scalar pV(0);
    scalar sV(0);
    forAll(primaryCells_, i)
    {
        const label& celli(primaryCells_[i]);
        pCOV += primV[celli]*primC[celli];
        pV += primV[celli];
    }
    reduce(pCOV, sumOp<vector>());
    reduce(pV, sumOp<scalar>());
    pCOV /= pV;
    forAll(secondaryCells_, i)
    {
        const label& celli(secondaryCells_[i]);
        sCOV += secV[celli]*secC[celli];
        sV += secV[celli];
    }
    reduce(sCOV, sumOp<vector>());
    reduce(sV, sumOp<scalar>());
    sCOV /= sV;
    vector delta(pCOV-sCOV);
    scalar errV((pV-sV)/pV);
    if ((pV-sV)/pV >= 1e-2)
    {
        FatalErrorInFunction
            << "Primary and secondary zone volumes differ by " << errV
            << exit(FatalError);
    }
    Info<< "Heat exchanger " << thisDictionary_.dictName() << " mapping: applied "
        << "translation of " << delta << " m" << endl;

    // Read mesh geometric data
    const pointField& pointsRef(sMesh_.points());
    const faceList& facesRef(sMesh_.faces());
    const cellList& cellsRef(sMesh_.cells());

    // Init new mesh geometric data
    pointField points(pointsRef.size());
    faceList faces(facesRef.size());
    cellList cells(cellsRef.size());

    // Create translated mesh
    // First, copy all points and apply translation
    forAll(points, i)
    {
        points[i] = pointsRef[i] + delta;
    }

    // Set faces, cells as copies of the original indexing (yet the indexing
    // applies to the translated points, so all the resulting faces and cells
    // will result translated as well)
    forAll(faces, i)
    {
        faces[i] = facesRef[i];
    }
    forAll(cells, i)
    {
        cells[i] = cellsRef[i];
    }

    // Assemble new translated mesh
    autoPtr<fvMesh> translatedMeshPtr;
    translatedMeshPtr.reset
    (
        new Foam::fvMesh
        (
            IOobject
            (
                "translatedMesh",
                mesh_.time().constant(),
                mesh_.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            std::move(points),
            std::move(faces),
            std::move(cells)
        )
    );
    fvMesh& translatedMesh = translatedMeshPtr();

    // Assemble mapping
    mappingPtr_.reset
    (
        new meshToMesh
        (
            translatedMesh,
            mesh_,
            Foam::meshToMesh::interpolationMethod
            (
                this->lookupOrDefault("interpolationMethod", 2)
            ),
            Foam::meshToMesh::procMapMethod::pmAABB,
            false
        )
    );
    // I do not care about the translatedMesh anymore after this point, I only
    // needed the mapping!
    translatedMeshPtr.clear();

    constructed = true;
}


// ************************************************************************* //
