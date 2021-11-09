/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
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

#include "massFlow.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(massFlow, 0);
    addToRunTimeSelectionTable(functionObject, massFlow, dictionary);
}
}

const Foam::Enum
<
    Foam::functionObjects::massFlow::regionType
>
Foam::functionObjects::massFlow::regionTypeNames_
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


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::massFlow::writeFileHeader(Ostream& os)
{
    if (writtenHeader_)
    {
        writeBreak(os);
    }
    else
    {
        writeHeader(os, "Field extents");
    }

    writeCommented(os, "Time");

    /*
    for (const word& fieldName : fieldSet_.selectionNames())
    {
        if (internalField_)
        {
            writeTabbed(os, fieldName + "_internal");
        }
        for (const label patchi : patchIDs_)
        {
            const word& patchName = mesh_.boundaryMesh()[patchi].name();
            writeTabbed(os, fieldName + "_" + patchName);
        }
    }
    */

    os  << endl;

    writtenHeader_ = true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::massFlow::massFlow
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict),
    regionName_(""),
    patchID_(0),
    faces_(0),
    alphaRhoPhiPtr_(nullptr),
    scaleFactor_(1.0)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::massFlow::read(const dictionary& dict)
{
    if (fvMeshFunctionObject::read(dict) && writeFile::read(dict))
    {
        regionType_ = regionType
        (
            regionTypeNames_.get
            (
                dict.get<word>
                (
                    "regionType"
                )
            )
        );
        regionName_ = dict.get<word>("regionName");
        alphaRhoPhiName_ = dict.get<word>("alphaRhoPhiName");
        scaleFactor_ = dict.lookupOrDefault<scalar>("scaleFactor", 1.0);

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

        return true;
    }

    return false;
}


bool Foam::functionObjects::massFlow::execute()
{
    return true;
}


bool Foam::functionObjects::massFlow::write()
{
    writeFileHeader(file());

    Log << type() << " " << name() <<  " write:" << nl;

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

    Log << "    " << regionTypeNames_[regionType_] << " " << regionName_ 
        << " massFlow = " << mDot << " kg/s over " << S << " m2" << endl;
    file() << mDot;
    this->setResult(regionName_+"_massFlow", mDot);

    Log << endl;

    return true;
}


// ************************************************************************* //
