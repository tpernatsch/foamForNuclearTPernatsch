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

#include "sliceMapperByMaterial.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * Static Data Member  * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sliceMapperByMaterial, 0);

    addToRunTimeSelectionTable
    (
        sliceMapper, 
        sliceMapperByMaterial, 
        dictionary
    );
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::sliceMapperByMaterial::calcAddressing() const
{
    // Initialize number of slices
    labelList matSliceCount(mat_.materialsList().size());
    labelList matSliceStart(mat_.materialsList().size(), 0);
    
    //- Read number of slices in each material
    nSlices_ = 0;
    
    forAll(mat_.materialsList(), i)
    {
        const materialModel& model = mat_.materialsList()[i];
        const dictionary& materialDict(model.materialModelDict());
        const label curSliceCount = materialDict.lookupOrDefault("nSlices", 1);
        
        matSliceCount[i] = curSliceCount;
        matSliceStart[i] = nSlices_;
        nSlices_ += curSliceCount;
    }
    
    // Allocate slices
    sliceAddrList_.clear();
    sliceAddrList_.setSize(nSlices_);
    isFuel_ = boolList(nSlices_, false);
    
    // Reference to sliceField
    scalarField& sliceFieldRef = sliceID_->ref();
    
    label curSlice = 0;
    
    forAll(mat_.materialsList(), i)
    {
        const label matSlices = matSliceCount[i];
        
        const labelList& addr = mat_.matAddrList()[i];
        const materialModel& model = mat_.materialsList()[i];
        const word& materialName = model.name();
        const dictionary& materialDict = model.materialModelDict();
        
        const vectorField& points = mesh_.points();
        const labelListList& cellPoints = mesh_.cellPoints();
        const vectorField& C = mesh_.C().internalField();

        const label curSliceStart = matSliceStart[i];
        
        // Tag current set of slices if this is a fuel material
        if (isA<fuelMaterial>(mat_.materialsList()[i]))
        {
            for(int m = 0; m < matSlices; m++)
            {
                label sliceI = curSliceStart + m;
                
                isFuel_[sliceI] = true;
            }
        }
        
        // Calculate slice heights
        //- Find min and max axial coordinate (along pinDirection_) in 
        //- current material.            
        scalar hMin(GREAT);
        scalar hMax(-GREAT);
        
        forAll(addr, j)
        {
            const label cellI = addr[j];

            //- Reference to cell as list of points
            const labelList& curCellPoints = cellPoints[cellI];

            forAll(curCellPoints, k)
            {   
                scalar hPoint = (points[curCellPoints[k]]) & pinDirection_;
                hMin = min(hMin, hPoint);
                hMax = max(hMax, hPoint);
            }
        }

        reduce(hMin, minOp<scalar>());
        reduce(hMax, maxOp<scalar>());
        scalar materialHeight = (hMax - hMin);

        //- Set height of each slice in current material
        scalarField sliceHeights(matSlices, materialHeight/matSlices);
        
        if (materialDict.found("heightSlices"))
        {
            sliceHeights = scalarField(materialDict.lookup("heightSlices"));

            //- Check if sum of slice heights correspond to hMax-hMin
            //- If not give error
            scalar sumHeights(sum(sliceHeights));
    
            if(sliceHeights.size() < matSlices)
            {
                FatalErrorIn("sliceMapperByMaterial")
                << "The number of slice heights given in \"heightSlices\" "
                << "is smaller than \"nSlices\"." << nl
                << "Check \"heightSlices\" in dictionary."
                << exit(FatalError);
            }
            else if(sliceHeights.size() > matSlices)
            {
                FatalErrorIn("sliceMapperByMaterial")
            << "The number of slice heights given in \"heightSlices\" "
                << "is larger than \"nSlices\"." << nl
                << "Check \"heightSlices\" in dictionary."
                << exit(FatalError);
            }          
            else if(mag(materialHeight - sumHeights) > 1e-6)
            {
                FatalErrorIn("sliceMapperByMaterial")
                << "The sum of slice heights ("
                << sumHeights << "m) "
                << "for material \"" << materialName << "\" "
                << "differs from calculated material height (" 
                << materialHeight << "m)." << nl
                << "Check \"heightSlices\" in dictionary."
                << exit(FatalError);
            }
        }

        // Assign correct slice ID to each cell
        // It's a local list i.e. maps the local cell ID (0 to addr.size()-1) 
        // to the local slice ID (0 to matSlices-1)
        labelList sliceIdsLocal(addr.size(), -1);
        
        forAll(addr, j)
        {
            const label cellI = addr[j];
            
            scalar hCell = C[cellI]&pinDirection_;
            
            //- Cell distance from bottom of cellZone
            scalar deltaH = hCell - hMin;

            //- Prepare sliceI corresponding to current cell
            label sliceI = 0;
            scalar totalHeight(0);   

            //- SliceN label for present cell
            for(int m = 0; m < matSlices; m++)
            {
                totalHeight += sliceHeights[m];

                if(totalHeight >= deltaH + 1e-6)
                {
                    sliceIdsLocal[j] = m;
                    sliceFieldRef[cellI] = m + curSliceStart;
                    break;
                }
            }
        }
        
        // Construct local slice addressing arrays for this material
        labelListList matSliceAddrListLocal = invertOneToMany(matSlices, sliceIdsLocal);
        
        for(int m = 0; m < matSlices; m++)
        {
            const label sliceI = m + curSliceStart;

            const labelList& localList = matSliceAddrListLocal[m];
            labelList globalList(localList.size());
            forAll(localList, k)
            {
                const label j = localList[k];
                globalList[k] = addr[j];
            }

            sliceAddrList_.set(sliceI, new labelList(globalList));

            // Check for empty slices
            scalar sizeI(sliceAddrList_[sliceI].size());
            reduce(sizeI, sumOp<scalar>());

            if(sizeI <= 0)
            {
                FatalErrorIn("sliceMapperByMaterial")
                << "Found empty slice (i.e. a slice with no cells) when "
                << "creating slices for material " << materialName << "." << nl
                << "Check the keyword \"nSlices\" in the corresponding materials "
                << "subdictionary within the solverDict." << nl
                << "As a general rule, the number of slices should be lower than " 
                << "the number of axial divisions in the mesh."
                << exit(FatalError);
            }
        }
    }

    sliceID_->correctBoundaryConditions();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sliceMapperByMaterial::sliceMapperByMaterial
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
    (
        mesh_.lookupObject<globalOptions>("globalOptions")
    );

    pinDirection_ = globalOpt.pinDirection();

    calcAddressing();   
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //



// ************************************************************************* //
