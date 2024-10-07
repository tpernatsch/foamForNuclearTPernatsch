/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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

#include "cellZoneMultiMaterialInterface.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cellZoneMultiMaterialInterface, 0);
    addToRunTimeSelectionTable
    (
        multiMaterialInterface, 
        cellZoneMultiMaterialInterface, 
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellZoneMultiMaterialInterface::cellZoneMultiMaterialInterface
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    uniformMultiMaterialInterface(mesh, dict),
    interfaceWeight_(readScalar(dict.lookup("interfaceWeights")))
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellZoneMultiMaterialInterface::~cellZoneMultiMaterialInterface()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

//- Return a surfaceScalarField of interface weights
void Foam::cellZoneMultiMaterialInterface::correct()
{
    //- Correct the default internal faces
    uniformMultiMaterialInterface::correct();
    
    //- Find zone interfaces
    labelList zoneIDs(mesh_.nCells(), -1);
    
    forAll(mesh_.cellZones(), zoneI)
    {
        IndirectList<label>(zoneIDs, mesh_.cellZones()[zoneI]) = zoneI;
    }
    
    {
        //- Internal faces at cellZone interfaces
        scalarField& w_i = weights_.ref();
        const labelList& owner = mesh_.owner();
        const labelList& neighbour = mesh_.neighbour();
        
        forAll(w_i, faceI)
        {
            if (zoneIDs[owner[faceI]] != zoneIDs[neighbour[faceI]])
            {
                w_i[faceI] = interfaceWeight_;
            }
        }
    }
    
    //- Update coupled boundaries
    forAll(weights_.boundaryField(), patchI)
    {
        fvsPatchField<scalar>& pf = weights_.boundaryFieldRef()[patchI];
        
        if(pf.coupled())
        {
            const lduInterface& p = refCast<const lduInterface>(pf.patch());

#ifdef OPENFOAMFOUNDATION            
            const unallocLabelList& owner = p.faceCells();
#elif OPENFOAMESI
            const UList<label>& owner = p.faceCells();
#endif            
            
            p.initInternalFieldTransfer(Pstream::defaultCommsType, zoneIDs);
            labelField neighbourZoneID = p.internalFieldTransfer
            (
                Pstream::defaultCommsType,
                zoneIDs
            );
            
            forAll(pf, faceI)
            {
                if (zoneIDs[owner[faceI]] != neighbourZoneID[faceI])
                {
                    pf[faceI] = interfaceWeight_;
                }
            }
        }
    }
}


// ************************************************************************* //
