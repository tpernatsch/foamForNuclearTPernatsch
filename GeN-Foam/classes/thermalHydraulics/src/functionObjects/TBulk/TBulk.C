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

#include "TBulk.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(TBulk, 0);
    addToRunTimeSelectionTable(functionObject, TBulk, dictionary);
}
}

const Foam::Enum
<
    Foam::functionObjects::TBulk::regionType
>
Foam::functionObjects::TBulk::regionTypeNames_
(
    {
        { 
            regionType::patch, 
            "patch" 
        },
        { 
            regionType::faceSet_, 
            "faceSet" 
        },
        { 
            regionType::faceZone_, 
            "faceZone" 
        }
    }
);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::TBulk::writeFileHeader(Ostream& os)
{
    if (writtenHeader_)
    {
        //writeBreak(os);
        os << endl;
    }
    else
    {
        writeHeader(os, "TBulk for regions: ");
        os << "    ";
        forAll(regionNames_, i)
        {
            os << regionNames_[i];
            if (i != regionNames_.size()-1)
                os << ", ";
            else
                os << endl;
        }
    }
    
    word time("Time = "+mesh_.time().timeName());

    os << time << endl;

    writtenHeader_ = true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::TBulk::TBulk
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict),
    regionNames_(0),
    patchIDs_(),
    thermoPtr_(nullptr),
    alphaRhoPhiPtr_(nullptr)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::TBulk::read(const dictionary& dict)
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
        regionNames_ = dict.get<wordReList>("regions");
        thermoName_ = dict.get<word>("thermo");
        alphaRhoPhiName_ = dict.get<word>("alphaRhoPhi");

        if (regionType_ == regionType::patch)
        {
            patchIDs_.clear();
            const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
            for (const wordRe& name : regionNames_)
            {
                patchIDs_.insert(pbm.findIndices(name));
            }
        }
        else if (regionType_ == regionType::faceSet_)
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
            for (const wordRe& name : regionNames_)
            {
                Info << name << " " << faceSets.found(name) << endl;
                if (faceSets.found(name))
                {
                    faceSet set(*faceSets[name]);
                    forAllIter(faceSet, set, iter)
                    {
                        label facei(*iter);
                        faces_.append(facei);
                    }
                }
            }
            
        }
        else 
        {
            FatalErrorInFunction
                << "faceZone region types currently not implemented"
                << exit(FatalError);
        }

        return true;
    }

    return false;
}


bool Foam::functionObjects::TBulk::execute()
{
    return true;
}


bool Foam::functionObjects::TBulk::write()
{
    writeFileHeader(file());

    Log << type() << " " << name() <<  " write:" << nl;

    //- Set pointers
    if (thermoPtr_ == nullptr)
    {
        thermoPtr_ = &mesh_.lookupObject<rhoThermo>(thermoName_);
    }
    if (alphaRhoPhiPtr_ == nullptr)
    {
        alphaRhoPhiPtr_ = &mesh_.lookupObject<surfaceScalarField>(alphaRhoPhiName_);
    } 
    
    //- Get refs for convenience
    const rhoThermo& thermo(*thermoPtr_);
    const volScalarField& T(thermo.T());
    const surfaceScalarField& alphaRhoPhi(*alphaRhoPhiPtr_);

    //- Cp is the only one that is computed on the fly, no refs to it
    tmp<volScalarField> Cp(thermo.Cp());

    switch (regionType_)
    {
        case regionType::patch :
        {   
            //- Total enthalpy flow (J/s) through patch
            scalar hDot(0.0);

            //- Total enthalpy flow per unit fluid temperature (J/s/K)
            scalar hDotByT(0.0);    
            
            for (const label patchi : patchIDs_)
            {
                const fvPatchScalarField& Cpp = Cp().boundaryField()[patchi];
                const fvPatchScalarField& Tp = T.boundaryField()[patchi];
                const fvsPatchField<scalar>& alphaRhoPhip = alphaRhoPhi.boundaryField()[patchi];
                forAll(Tp, i)
                {
                    scalar alphaRhoCpUDotSf
                    (
                        //alphap[i]*rhop[i]*Cpp[i]*Up[i]&Sf[i]
                        alphaRhoPhip[i]*Cpp[i]
                    );
                    hDotByT += alphaRhoCpUDotSf;
                    hDot += alphaRhoCpUDotSf*Tp[i];
                }
            }
            reduce(hDot, sumOp<scalar>());
            reduce(hDotByT, sumOp<scalar>());
            scalar TBulkValue(hDot/max(hDotByT, 1e-6));

            Log << "    patch(es) " << regionNames_ << " TBulk = " << TBulkValue
                << " K" << endl;
            file() << "TBulk = " << TBulkValue << " K";
            //this->setResult(compoundPatchesName_ + "_TBulk", TBulkValue);
            break;
        }
        case regionType::faceSet_ :
        {
            //- Total enthalpy flow (J/s) through faceSet
            scalar hDot(0.0);

            //- Total enthalpy flow per unit fluid temperature (J/s/K)
            scalar hDotByT(0.0); 
            surfaceScalarField Cpp(fvc::interpolate(Cp));
            surfaceScalarField Tp(fvc::interpolate(T));

            forAll(faces_, i)
            {
                const label& facei(faces_[i]);
                scalar alphaRhoCpUDotSf
                (
                    alphaRhoPhi[facei]*Cpp[facei]
                );
                hDotByT += alphaRhoCpUDotSf;
                hDot += alphaRhoCpUDotSf*Tp[facei];
            }

            reduce(hDot, sumOp<scalar>());
            reduce(hDotByT, sumOp<scalar>());
            scalar TBulkValue(hDot/max(hDotByT, 1e-6));

            Log << "    faceSet(s) " << regionNames_ << " TBulk = " << TBulkValue
                << " K" << endl;
            file() << "TBulk = " << TBulkValue << " K";

            break;
        }
        default :
        {
            FatalErrorInFunction
                << "faceZone region types currently not implemented"
                << exit(FatalError);
            break;
        }
    }

    Log << endl;

    return true;
}


// ************************************************************************* //
