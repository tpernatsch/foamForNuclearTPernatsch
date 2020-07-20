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

#include "massFlowFlux.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(massFlowFlux, 0);
    addToRunTimeSelectionTable(functionObject, massFlowFlux, dictionary);
}
}

const Foam::Enum
<
    Foam::functionObjects::massFlowFlux::regionType
>
Foam::functionObjects::massFlowFlux::regionTypeNames_
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

void Foam::functionObjects::massFlowFlux::writeFileHeader(Ostream& os)
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

Foam::functionObjects::massFlowFlux::massFlowFlux
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict),
    patchNames_(0),
    patchIDs_(),
    compoundPatchesName_(""),
    thermoPtr_(nullptr),
    alphaPtr_(nullptr),
    UPtr_(nullptr)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::massFlowFlux::read(const dictionary& dict)
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
        thermoName_ = dict.get<word>("thermo");
        alphaName_ = dict.lookupOrDefault<word>("alpha", "none");

        if (regionType_ == regionType::patch)
        {
            patchIDs_.clear();
            const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
            patchNames_ = dict.get<wordReList>("patches");
            for (const wordRe& name : patchNames_)
            {
                patchIDs_.insert(pbm.findIndices(name));
                compoundPatchesName_ += name+":";
            }
            UName_ = dict.get<word>("U");
        }
        else
        {
            FatalErrorInFunction
                << "faceSet, faceZone region types currently not implemented"
                << exit(FatalError);
        }

        return true;
    }

    return false;
}


bool Foam::functionObjects::massFlowFlux::execute()
{
    return true;
}


bool Foam::functionObjects::massFlowFlux::write()
{
    writeFileHeader(file());

    Log << type() << " " << name() <<  " write:" << nl;

    //- Set pointers
    if (thermoPtr_ == nullptr)
    {
        thermoPtr_ = &mesh_.lookupObject<rhoThermo>(thermoName_);
    }
    if (UPtr_ == nullptr)
    {
        UPtr_ = &mesh_.lookupObject<volVectorField>(UName_);
    }
    if (alphaPtr_ == nullptr)
    {
        if (alphaName_ != "none")
        {
            alphaPtr_ = &mesh_.lookupObject<volScalarField>(alphaName_);
        }
        else
        {
            alphaPtr_ = new volScalarField
            (
                IOobject
                (
                    "oneField",
                    mesh_.time().timeName(),
                    mesh_
                ),
                mesh_,
                dimensionedScalar("", dimless, 1.0)
            );
        }
    }
    
    //- Get refs for convenience
    const rhoThermo& thermo(*thermoPtr_);
    const volVectorField& U(*UPtr_);
    const volScalarField& alpha(*alphaPtr_);
    const volScalarField& rho(thermo.rho());

    switch (regionType_)
    {
        case regionType::patch :
        {   
            //- Total patch surface
            scalar S(0.0); 

            //- Total mass flow through patch
            scalar VDot(0.0);

            //- Total mass flux through patch
            scalar mDot(0.0);
            
            for (const label patchi : patchIDs_)
            {
                const fvPatchVectorField& Up = U.boundaryField()[patchi];
                const fvPatchScalarField& alphap = 
                    alpha.boundaryField()[patchi];
                const fvPatchScalarField& rhop = rho.boundaryField()[patchi];
                const fvPatch& patch(mesh_.boundary()[patchi]);
                const scalarField& magSf(patch.magSf());
                const vectorField& Sf(patch.Sf());
                forAll(Sf, i)
                {
                    const scalar& magSfi(magSf[i]);
                    scalar alphaUDotSfi(alphap[i]*Up[i]&Sf[i]);
                    S += magSfi;
                    VDot += alphaUDotSfi;
                    mDot += rhop[i]*alphaUDotSfi;
                }
            }
            reduce(S, sumOp<scalar>());
            reduce(VDot, sumOp<scalar>());
            reduce(mDot, sumOp<scalar>());

            scalar mFlux(mDot/S);

            Log << "    patch(es) " << patchNames_ << " tot Surf = " << S 
                << " m2" << nl
                << " Flow = " << VDot << " m3/s " << nl
                << " massFlow = " << mDot << " kg/s " << nl
                << " massFlux = " << mFlux << " kg/s/m2" << endl;
            file() << mFlux;
            this->setResult(compoundPatchesName_ + "_massFlowFlux", mFlux);

            break;
        }
        default :
        {
            FatalErrorInFunction
                << "faceSet, faceZone region types currently not implemented"
                << exit(FatalError);
            break;
        }
    }

    Log << endl;

    return true;
}


// ************************************************************************* //
