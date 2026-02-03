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

#include "sliceMapperByPellets.H"
#include "addToRunTimeSelectionTable.H"
#include "SortableList.H"


// * * * * * * * * * * * * * Static Data Member  * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sliceMapperByPellets, 0);

    addToRunTimeSelectionTable
    (
        sliceMapper, 
        sliceMapperByPellets, 
        dictionary
    );
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::sliceMapperByPellets::calcAddressing() const
{
    // Initialize number of slices
    nSlices_ = 0;

    sliceAddrList_.clear();
    const labelListList& matAddrList(mat_.matAddrList());

    // Reference to sliceField
    scalarField& sliceFieldRef = sliceID_->ref();

    forAll(mat_.materialsList(), i)
    {
        const labelList& addrI(matAddrList[i]);

        dictionary materialDict(mat_.materialsList()[i].materialModelDict());
        const word materialName(mat_.materialsList()[i].name());        

        // Find min and max axial coordinate (along pinDirection_) in 
        // current material.            
        scalar hMin(GREAT);
        scalar hMax(-GREAT);
                
        forAll(addrI, j)
        {
            const label cellI = addrI[j];

            // Reference to cell as list of points
            const labelList& c = mesh_.cellPoints()[cellI]; 

            forAll(c, pointI)
            {   
                // Reference to point
                const vector& p = mesh_.points()[c[pointI]];
                
                scalar hPoint(p&pinDirection_);
                hMin = min(hMin, hPoint);
                hMax = max(hMax, hPoint);
            }
        }

        reduce(hMin, minOp<scalar>());
        reduce(hMax, maxOp<scalar>());
        scalar materialHeight = (hMax - hMin);

        // Read number of pellets in current material (1 for non-fuel materials)
        scalar nPellets(isA<fuelMaterial>(mat_.materialsList()[i])?
            readInt(materialDict.lookup("nPellets")) : 1);

        // Set height of each pellet in current material
        scalarField sliceHeights(nPellets, materialHeight/nPellets);

        // Reference to cell centers
        const volVectorField& C = mesh_.C();

        // Slice addressing for current material
        List<labelList> sliceAddrListI(nPellets, List<label>(0, 0));

        // Assign correct slice ID to each cell and collect all cell 
        // addresses for each slice in cellIDsPerSlice_
        forAll(addrI, j)
        {
            const label cellI = addrI[j];
            
            scalar hCell = C[cellI]&pinDirection_;
            
            // Cell distance from bottom of cellZone
            scalar deltaH = hCell - hMin;

            // Prepare sliceID corresponding to current cell
            scalar currentSliceID(0);   
            scalar totalHeight(0);   

            // SliceN label for present cell
            for(int m=0; m<nPellets; m++)
            {
                totalHeight += sliceHeights[m];

                if(totalHeight >= deltaH + 1e-6)
                {
                    currentSliceID = m;
                    break;
                }
            }
            
            sliceAddrListI[currentSliceID].append(cellI);
        }         

        // If the cells are taller than the pellet size, some of the slices
        // will be empty. The following lines pass only the non-empty slices
        // to the global sliceAddrList_
        label matSlices(0);
        forAll(sliceAddrListI, matSliceID)
        {
            scalar sizeI(sliceAddrListI[matSliceID].size());
            reduce(sizeI, sumOp<scalar>());

            // Create the slice only if the size is bigger than 0 in at least
            // 1 processor (the local ID list might be empty)
            if(sizeI>0)
            {
                sliceAddrList_.append(new labelList(sliceAddrListI[matSliceID]));
                isFuel_.append(isA<fuelMaterial>(mat_.materialsList()[i]));
                forAll(sliceAddrList_[nSlices_ + matSlices], i)
                {
                    label cellI = sliceAddrList_[nSlices_ + matSlices][i];
                    sliceFieldRef[cellI] = nSlices_ + matSlices;
                }
                matSlices += 1;
            }
        }

        // Check for empty slices
        forAll(sliceAddrList_, sliceI)
        {
            scalar sizeI(sliceAddrList_[sliceI].size());
            reduce(sizeI, sumOp<scalar>());

            if(sizeI<=0)
            {
                FatalErrorIn("sliceMapperByPellets")
                << "Found empty slice (i.e. a slice with no cells) when "
                << "creating slices for material " << materialName << "." << nl
                << "Check the keyword \"nSlices\" in the corresponding materials "
                << "subdictionary within the solverDict." << nl
                << "As a general rule, the number of slices should be lower than " 
                << "the number of axial divisions in the mesh."
                << exit(FatalError);
            }
        }

        // Update total number of slice
        nSlices_ += matSlices;      
    }

    sliceID_->correctBoundaryConditions();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sliceMapperByPellets::sliceMapperByPellets
(
    const fvMesh& mesh,
    const materials& mat,
    const dictionary& sliceMapperDict
)
:
    sliceMapper(mesh, mat, sliceMapperDict)
{   
    createSliceID();
    
    const globalOptions& globalOpt
    (mesh_.lookupObject<globalOptions>("globalOptions"));

    pinDirection_ = (globalOpt.pinDirection());

    calcAddressing();
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// ************************************************************************* //
