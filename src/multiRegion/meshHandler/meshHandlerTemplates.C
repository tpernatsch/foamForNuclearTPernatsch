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

#include "meshHandler.H"
#include "solver.H"
#include "Time.H"
#include "interpolationCellPoint.H"
#include "radialBasisFunctionInterpolation.H"
#include "mergeOrSplitBaffles.H"
#include "hexCellFvMesh.H"
#include "dynamicFvMesh.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class Type>
void Foam::meshHandler::map
(
    word srcFieldName,
    word tgtFieldName,
    label tgtFieldLabel,
    label srcFieldLabel,
    bool removeBaffles
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    if(this->meshes_[tgtFieldLabel].foundObject<VolFieldType>(tgtFieldName))
    {

        if
        (
            (
                removeBaffles
                && !(this->mappingMeshes_[srcFieldLabel].foundObject<VolFieldType>(srcFieldName))
            ) || (
                !removeBaffles
                && !(this->meshes_[srcFieldLabel].foundObject<VolFieldType>(srcFieldName))
            )
        )
        {
            FatalErrorInFunction
            << "Types of mapping fields " << srcFieldName <<
            " and " << tgtFieldName << " do not match" <<endl;
        }

        VolFieldType& tgtField = const_cast<VolFieldType&>
        (
            this->meshes_[tgtFieldLabel].lookupObject<VolFieldType>(tgtFieldName)
        );


        mappingList_[tgtFieldLabel][srcFieldLabel].mapTgtToSrc
        (
            removeBaffles
                ? this->mappingMeshes_[srcFieldLabel].lookupObject<VolFieldType>(srcFieldName)
                : this->meshes_[srcFieldLabel].lookupObject<VolFieldType>(srcFieldName),
            plusEqOp<Type>(),
            tgtField
        );
        tgtField.correctBoundaryConditions();
    }
}

template<class Type>
void Foam::meshHandler::mapAndWrite
(
    word srcFieldName,
    word tgtFieldName,
    label tgtFieldLabel,
    label srcFieldLabel,
    bool removeBaffles,
    bool initializeMappedFields
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    if(this->meshes_[tgtFieldLabel].foundObject<VolFieldType>(tgtFieldName))
    {

        if
        (
            (
                removeBaffles
                && !(this->mappingMeshes_[srcFieldLabel].foundObject<VolFieldType>(srcFieldName))
            ) || (
                !removeBaffles
                && !(this->meshes_[srcFieldLabel].foundObject<VolFieldType>(srcFieldName))
            )
        )
        {
            FatalErrorInFunction
            << "Types of mapping fields " << srcFieldName <<
            " and " << tgtFieldName << " do not match" <<endl;
        }

        VolFieldType& tgtField = const_cast<VolFieldType&>
        (
            this->meshes_[tgtFieldLabel].lookupObject<VolFieldType>(tgtFieldName)
        );

        mappingList_[tgtFieldLabel][srcFieldLabel].mapTgtToSrc
        (
            removeBaffles
                ? this->mappingMeshes_[srcFieldLabel].lookupObject<VolFieldType>(srcFieldName)
                : this->meshes_[srcFieldLabel].lookupObject<VolFieldType>(srcFieldName),
            plusEqOp<Type>(),
            tgtField
        );
        tgtField.correctBoundaryConditions();
        if (initializeMappedFields)
        {
            tgtField.write();
        }
    }
}


template<class Type>
void Foam::meshHandler::createBaffleLessField(word fieldName, label regionLabel, fvMesh& baffleLessMesh, label& nCouplingFields)
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    // First: check if the field exist on the normal mesh for template reasons

    if(meshes_[regionLabel].foundObject<VolFieldType>(fieldName))
    {
        // Check that field has not been created yet
        if(!(baffleLessMesh.foundObject<VolFieldType>(fieldName+".baffleLess")))
        {
            // Create the cpupling field
            couplingFields_[regionLabel].resize(nCouplingFields+1);
            couplingFields_[regionLabel][nCouplingFields]=
            new VolFieldType
            (
                IOobject
                (
                    IOobject::groupName(fieldName, "baffleLess"),
                    baffleLessMesh.time().name(),
                    baffleLessMesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                baffleLessMesh,
                dimensioned<Type>("", meshes_[regionLabel].lookupObject<VolFieldType>(fieldName).dimensions(), Zero),
                zeroGradientFvPatchScalarField::typeName
            );
            nCouplingFields++;

            Info << "Creating baffleless field "<< fieldName<<".baffleLess" <<endl;
        }
    }

}

// ************************************************************************* //
