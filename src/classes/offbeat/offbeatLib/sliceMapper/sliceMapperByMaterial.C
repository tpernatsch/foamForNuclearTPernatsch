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

        //- Find min and max axial coordinate (along pinDirection_) in 
        //- current material.            
        scalar hMin(GREAT);
        scalar hMax(-GREAT);
                
        forAll(addrI, j)
        {
            const label cellI = addrI[j];

            //- Reference to cell as list of points
            const labelList& c = mesh_.cellPoints()[cellI]; 

            forAll(c, pointI)
            {   
                //- Reference to point
                const vector& p = mesh_.points()[c[pointI]];
                
                scalar hPoint(p&pinDirection_);
                hMin = min(hMin, hPoint);
                hMax = max(hMax, hPoint);
            }
        }

        reduce(hMin, minOp<scalar>());
        reduce(hMax, maxOp<scalar>());
        scalar materialHeight = (hMax - hMin);

        //- Read number of slices in current material
        scalar matSlices(readInt(materialDict.lookup("nSlices")));

        //- Set height of each slice in current material
        scalarField sliceHeights = 
        ( 
            materialDict.found("heightSlices")
        ) ?        
        materialDict.lookup("heightSlices")
        : 
        scalarField(matSlices, materialHeight/matSlices);

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

        //- Reference to cell centers
        const volVectorField& C = mesh_.C();

        //- Add an empty lists per slice to sliceAddrList_
        for(int j=0; j<matSlices; j++)
        {
            sliceAddrList_.append(new List<label>(0, 0));  
            isFuel_.append(isA<fuelMaterial>(mat_.materialsList()[i]));    
        }   

        //- Assign correct slice ID to each cell and collect all cell 
        //- addresses for each slice in cellIDsPerSlice_
        forAll(addrI, j)
        {
            const label cellI = addrI[j];
            
            scalar hCell = C[cellI]&pinDirection_;
            
            //- Cell distance from bottom of cellZone
            scalar deltaH = hCell - hMin;

            //- Prepare sliceID corresponding to current cell
            scalar sliceID(0);   
            scalar totalHeight(0);   

            //- SliceN label for present cell
            for(int m=0; m<matSlices; m++)
            {
                totalHeight += sliceHeights[m];

                if(totalHeight >= deltaH + 1e-6)
                {
                    //- Add also the total number of slice
                    sliceID = m + nSlices_;
                    break;
                }
            }
            
            sliceFieldRef[cellI] = sliceID;
            sliceAddrList_[sliceID].append(cellI);

        }         

        // Check for empty slices
        forAll(sliceAddrList_, sliceI)
        {
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

        //- Update total number of slice
        nSlices_ += matSlices;        
    }

    sliceID_->correctBoundaryConditions();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sliceMapperByMaterial::sliceMapperByMaterial
(
    const fvMesh& mesh,
    const materials& mat,
    const dictionary& mapperOptDict
)
:
    sliceMapper(mesh, mat, mapperOptDict)
{   
    createSliceID();
    
    const globalOptions& globalOpt
    (mesh_.lookupObject<globalOptions>("globalOptions"));

    pinDirection_ = (globalOpt.pinDirection());

    calcAddressing();   
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //



// ************************************************************************* //
