/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2412                                                  |
|    Copyright 2011-2016 OpenFOAM Foundation, 2017-2024 OpenCFD Ltd.          |
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
#if defined __has_include
#  if __has_include(<commDataLayer.H>)
#    include <commDataLayer.H>
#    define isCommDataLayerIncluded
#  endif
#endif

#ifdef isCommDataLayerIncluded

#include "massFlowToFMU.H"
#include "addToRunTimeSelectionTable.H"
#include "commDataLayer.H"
#include "externalIOObject.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
/*
// public fvMeshFunctionObject  //This would also work
// but one would have to put the functionObject in
// the controlDict instead of the externalCouplingDict
// which goes agaist the logif of FMU4FOAM
namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(massFlowToFMU, 0);
    addToRunTimeSelectionTable(functionObject, massFlowToFMU, dictionary);
}
}
*/

namespace Foam
{
namespace externalIOObject
{
    defineTypeNameAndDebug(massFlowToFMU, 0);
    addToRunTimeSelectionTable(externalIOObject, massFlowToFMU, dictionary);
}
}

const Foam::Enum
<
    Foam::externalIOObject::massFlowToFMU::regionType
>
Foam::externalIOObject::massFlowToFMU::regionTypeNames_
(
    {
        {
            regionType::patch,
            "patch"
        },
        {
            regionType::faceSet,
            "faceSet"
        },
        {
            regionType::faceZone,
            "faceZone"
        }
    }
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::externalIOObject::massFlowToFMU::massFlowToFMU
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    externalIOObject(name, runTime, dict),
    mesh_
    (
        refCast<const fvMesh>
        (
            time_.lookupObject<objectRegistry>
            (
                dict.getOrDefault("region", polyMesh::defaultRegion)
            )
        )
    ),
    regionName_(dict.get<word>("regionName")),
    patchID_(0),
    faces_(0),
    alphaRhoPhiName_(dict.get<word>("alphaRhoPhiName")),
    alphaRhoPhiPtr_(nullptr),
    scaleFactor_(dict.lookupOrDefault<scalar>("scaleFactor", 1.0)),
    nameFMU_(dict.get<word>("nameFMU")),
    initValue_(dict.lookupOrDefault<scalar>("initValue", 0.0))
{
    read(dict);
    execute();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool Foam::externalIOObject::massFlowToFMU::read(const dictionary& dict)
{
    // Set region type
    regionType_ = regionType
    (
        regionTypeNames_.get
        (
            dict.get<word>("regionType")
        )
    );

    if (regionType_ == regionType::patch)
    {
        const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
        patchID_ = pbm.findPatchID(regionName_);
    }
    else if (regionType_ == regionType::faceZone)
    {
        const faceZoneMesh& faceZones(mesh_.faceZones());
        const labelList& faceis(faceZones[regionName_]);
        forAll(faceis, i)
        {
            faces_.append(faceis[i]);
        }
    }
    else if (regionType_ == regionType::faceSet)
    {
        IOobjectList objects
        (
            mesh_,
            mesh_.time().findInstance
            (
                polyMesh::meshSubDir/"sets",
                word::null,
                IOobject::READ_IF_PRESENT,
                mesh_.facesInstance()
            ),
            polyMesh::meshSubDir/"sets"
        );
        IOobjectList faceSets(objects.lookupClass(faceSet::typeName));
        if (faceSets.found(regionName_))
        {
            Foam::faceSet set(*faceSets[regionName_]);
            forAllIter(faceSet, set, iter)
            {
                label facei(*iter);
                faces_.append(facei);
            }
        }
    }

    // Set FMI initialization
    commDataLayer& data = commDataLayer::New(time_);

    data.storeObj(initValue_, nameFMU_, commDataLayer::causality::out);

    return false;
}

bool Foam::externalIOObject::massFlowToFMU::execute()
{
    commDataLayer& data = commDataLayer::New(time_);

    scalar& result = data.getObj<scalar>(
        nameFMU_,
        commDataLayer::causality::out
    );

    if (alphaRhoPhiPtr_ == nullptr)
    {
        alphaRhoPhiPtr_ =
            &mesh_.lookupObject<surfaceScalarField>(alphaRhoPhiName_);
    }

    const surfaceScalarField& alphaRhoPhi(*alphaRhoPhiPtr_);

    scalar S(0.0);
    scalar mDot(0.0);

    if (regionType_ == regionType::patch)
    {
        const fvsPatchField<scalar>& alphaRhoPhip
            = alphaRhoPhi.boundaryField()[patchID_];
        const fvPatch& patch(mesh_.boundary()[patchID_]);
        const scalarField& magSf(patch.magSf());
        forAll(magSf, i)
        {
            const scalar& magSfi(magSf[i]);
            S += magSfi;
            mDot += mag(alphaRhoPhip[i]);
        }
    }
    else
    {
        const scalarField& magSf(mesh_.magSf());
        forAll(faces_, i)
        {
            const label& facei(faces_[i]);
            const scalar& magSfi(magSf[facei]);
            S += magSfi;
            mDot += mag(alphaRhoPhi[facei]);
        }
    }

    reduce(S, sumOp<scalar>());
    reduce(mDot, sumOp<scalar>());

    mDot *= scaleFactor_;

    result = mDot;

    return false;
}

bool Foam::externalIOObject::massFlowToFMU::write()
{
    return false;
}

#endif

// ************************************************************************* //
