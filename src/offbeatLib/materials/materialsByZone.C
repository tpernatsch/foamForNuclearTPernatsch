/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "materialsByZone.H"
#include "zeroGradientFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "emptyFvPatchFields.H"
#include "fuelMaterial.H"
#include "processorFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(materialsByZone, 0);
    addToRunTimeSelectionTable
    (
        materials, 
        materialsByZone, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::materialsByZone::checkZoneDefinition()
{
    labelList zoneMap(mesh_.nCells(), -1);
    
    forAll(mesh_.cellZones(), zoneI)
    {
        const labelList& addr = mesh_.cellZones()[zoneI];
        
        forAll(addr, i)
        {
            const label cellI = addr[i];
            
            if (zoneMap[cellI] != -1)
            {
                FatalErrorIn("Foam::materialsByZone::checkZoneDefinition()")
                    << "Cell zones " << mesh_.cellZones()[zoneI].name()
                    << " and " << mesh_.cellZones()[zoneMap[cellI]].name()
                    << " overlap" << endl << Foam::abort(FatalError);
            }
            zoneMap[cellI] = zoneI;
        }
    }
    
    forAll(zoneMap, cellI)
    {
        if (zoneMap[cellI] == -1)
        {
            FatalErrorIn("Foam::materialsByZone::checkZoneDefinition()")
                << "Cell zones do not include all mesh cells." << endl
                << Foam::abort(FatalError);
        }
    }
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::materialsByZone::calcAddressing(const dictionary& materialsDict) const
{
    matAddrList_.setSize(mesh_.cellZones().size());
    forAll(matAddrList_, i)
    {
        matAddrList_[i] = mesh_.cellZones()[i];
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::materialsByZone::materialsByZone
(
    const fvMesh& mesh, 
    const dictionary& materialsDict
)
:
    materials(mesh, materialsDict)
{
    matList_.clear();
    matList_.resize(mesh_.cellZones().size());

    checkZoneDefinition();

    calcAddressing(materialsDict);

    forAll(matList_, zoneI)
    {
        const word zoneName = mesh_.cellZones()[zoneI].name();

        const dictionary& materialModelDict(materialsDict.subDict(zoneName));
        
        Info << "     For cellzone " << zoneName << ": " << endl;

        matList_.set
        (
            zoneI,
            materialModel::New
            (
                mesh_,  
                materialModelDict,
                matAddrList_[zoneI]
            )
        );

        Info << endl;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::materialsByZone::~materialsByZone()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::materialsByZone::correctThermoMechProperties()
{    
    if(T_ == nullptr)
    {
        T_ = &mesh_.lookupObject<volScalarField>("T");
    }

    // Correct material properties
    const scalarField& Ti = T_->internalField();
    scalarField& rhoi = rho_.ref();
    scalarField& Cpi = Cp_.ref();
    scalarField& ki = k_.ref();
    scalarField& epsi = emissivity_.ref();
    scalarField& Ei = E_.ref();
    scalarField& nui = nu_.ref();

    forAll(matList_, i)
    {
        materialModel& matI(matList_[i]);

        matI.rho(rhoi, Ti);
        matI.Cp(Cpi, Ti);
        matI.k(ki, Ti);
        matI.emissivity(epsi, Ti);
        matI.E(Ei, Ti);
        matI.nu(nui, Ti);
        matI.soften(Ei, nui);
    }

    updateDensity();

    // Lame's parameters
    lambda_ = nu_*E_/((1.0 + nu_)*(1.0 - 2.0*nu_));
    mu_ = E_/(2.0*(1.0 + nu_));
    threeK_ = E_/(1.0 - 2.0*nu_);
    
    //- Needed if we set boundary conditions for k and emmissivity to calculated
    forAll(mesh_.boundaryMesh(), patchI)
    {
        fvPatchScalarField& kp = k_.boundaryFieldRef()[patchI];
        kp.operator=(kp.patchInternalField());

        fvPatchScalarField& emmisivity_p = emissivity_.boundaryFieldRef()[patchI];
        emmisivity_p.operator=(emmisivity_p.patchInternalField());
    }

    rho_.correctBoundaryConditions();
    Cp_.correctBoundaryConditions();
    k_.correctBoundaryConditions();
    emissivity_.correctBoundaryConditions();
    E_.correctBoundaryConditions();
    nu_.correctBoundaryConditions();
    lambda_.correctBoundaryConditions();
    mu_.correctBoundaryConditions();
    threeK_.correctBoundaryConditions();
}


void Foam::materialsByZone::correctBehavioralModels()
{  
    if(T_ == nullptr)
    {
        T_ = &mesh_.lookupObject<volScalarField>("T");
    }

    const scalarField& Ti = T_->internalField();  
    
    forAll(matList_, i)
    {
        materialModel& matI(matList_[i]);
        matI.correctBehavioralModels(Ti);
    }
}


void Foam::materialsByZone::checkFailure()
{  
    forAll(matList_, i)
    {
        materialModel& matI(matList_[i]);
        matI.checkFailure();
    }
}


Foam::tmp<Foam::volSymmTensorField>
Foam::materialsByZone::thermalExpansion() const
{    
    tmp<volSymmTensorField> teps
    (
        new volSymmTensorField
        (
            IOobject
            (
                "threeK*alpha*T",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedSymmTensor("threeK*alpha*T", dimless, symmTensor::zero)
        )
    );
    
    symmTensorField& epsi = teps->ref();
        
    const scalarField& Ti = T_->internalField();

    // We want to update the thermal expansion also on the patches (calculated
    // BC instead of zeroGradient). Because different materials can share the
    // same boundary (e.g. UO2 and UPuO2 in the same simulation sharing the 
    // same fuelOuter surface), the following loop is necessary. 
    // For each material, it finds the enclosing patches, and for each patch
    // finds the appropriate list of faceIDS.
    forAll(matList_, i)
    {
        const labelList& addr = mesh_.cellZones()[i];
        labelList patchIDs_(0);
        PtrList<List<label>> addrList(0);

        // Loop over all cells of current material
        forAll(addr, addrI)
        {
            const label cellI = addr[addrI];

            // Reference to list of faces per cell
            const cell& c = mesh_.cells()[cellI];     

            // Loop over all faces of current cell
            forAll(c, i)
            {
                const label patchID = mesh_.boundaryMesh().whichPatch(c[i]);
                label indexInPatchID(0);

                if (patchID > -1)
                {
                    bool found_ = false;

                    // Check wheter patchID has been added to list of patchIDs 
                    // for current material
                    forAll(patchIDs_, patchi)
                    {
                        indexInPatchID = patchi;
                        if (patchID ==  patchIDs_[patchi])
                        {
                            found_ = true;
                            break;
                        }
                    }

                    // Add patch and face IDs only for non-empty and non-proc
                    // patches
                    if
                    (
                        !isType<processorFvPatchScalarField>
                        (
                            T_->boundaryField()[patchID]
                        )
                        and                        
                        !isType<emptyFvPatchScalarField>
                        (
                            T_->boundaryField()[patchID]
                        )
                    )
                    {

                        // Add patchID to list for current material
                        if(!found_)
                        {
                            patchIDs_.append(patchID);
                            indexInPatchID = addrList.size();
                            addrList.append(new List<label>(0, 0)); 
                        }

                        // Add faceID to the list for the current patch for the 
                        // current material
                        label faceID = 
                        mesh_.boundaryMesh()[patchID].whichFace(c[i]);
                        addrList[indexInPatchID].append(faceID);
                    }

                }
            }
        }

        // NOTE: necessary for implicitContact?
        forAll(patchIDs_, patchI)
        {
            const label patchID = patchIDs_[patchI];
            const scalarField& TP = T_->boundaryField()[patchID];

            symmTensorField& epsP = teps->boundaryFieldRef()[patchID];
            matList_[i].alphaT(epsP, TP, addrList[patchI]);
        }

        matList_[i].alphaT(epsi, Ti, addr);
    }

    teps->correctBoundaryConditions();
    
    return teps;
}


Foam::tmp<Foam::volSymmTensorField>
Foam::materialsByZone::additionalStrains() const
{
    tmp<volSymmTensorField> tempStrain
    (
        new volSymmTensorField
        (
            IOobject
            (
                "additionalStrain",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedSymmTensor("additionalStrain", dimless, symmTensor::zero),
            "zeroGradient"
        )
    );
    
    symmTensorField& epsi = tempStrain->ref();        
    
    forAll(matList_, i)
    {
        const materialModel& matI(matList_[i]);
        matI.additionalStrains(epsi);
    }

    tempStrain->correctBoundaryConditions();
    
    return tempStrain;
}

// ************************************************************************* //
