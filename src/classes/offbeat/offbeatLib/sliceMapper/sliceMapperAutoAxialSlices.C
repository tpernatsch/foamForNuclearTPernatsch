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

#include "sliceMapperAutoAxialSlices.H"
#include "addToRunTimeSelectionTable.H"
#include "SortableList.H"


// * * * * * * * * * * * * * Static Data Member  * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sliceMapperAutoAxialSlices, 0);

    addToRunTimeSelectionTable
    (
        sliceMapper, 
        sliceMapperAutoAxialSlices, 
        dictionary
    );
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::sliceMapperAutoAxialSlices::calcAddressing() const
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

        // Cell centers axial position
        const scalarField Cz = mesh_.C()&pinDirection_;

        // HashTable for storing the slice addressing for current material
        HashTable<labelList, scalar> sliceAddrListI(1);

        // Create an entry of the hash table for each different axial location 
        // and append the cellID
        forAll(addrI, j)
        {
            const label cellI = addrI[j];
            
            scalar CzI = round(Cz[cellI]/precision_)*precision_;

            HashTable<labelList, scalar>::iterator iter = sliceAddrListI.find(CzI);

            if (iter != sliceAddrListI.end())
            {
                iter().append(cellI);
            }
            else
            {
                sliceAddrListI.insert(CzI, labelList(1, cellI));
            }
        }

        // Handle multiprocessor simulations //

        // Gather HashTable.toc for each processor in one list 
        // (one height list per processor) and then redistribute
        List< List<scalar> > gatheredAxialLocationLists(Pstream::nProcs());
        gatheredAxialLocationLists[Pstream::myProcNo()] = sliceAddrListI.toc();
        Pstream::gatherList(gatheredAxialLocationLists);
        Pstream::scatterList(gatheredAxialLocationLists);

        // Create a global list of axial location in each processor and sort it
        // (some locations might be repeated in different processors)
        std::set<scalar> sortedAxialLocationList;
        forAll(gatheredAxialLocationLists, procI)
        {
            scalarList procList(gatheredAxialLocationLists[procI]);
            forAll(procList, hI)
            {
                sortedAxialLocationList.insert(procList[hI]);
            }
        }        

        // Create slices for current material keeping track of global number of 
        // slices (summing across the various materials)
        label matSlices(0);
        for(const scalar& zI :  sortedAxialLocationList)
        {
            // Append list of IDs
            // NOTE: in multiprocessor the list will be empty
            sliceAddrList_.append(new labelList(sliceAddrListI(zI)));
            isFuel_.append(isA<fuelMaterial>(mat_.materialsList()[i]));
            forAll(sliceAddrList_[nSlices_ + matSlices], i)
            {
                label cellI = sliceAddrList_[nSlices_ + matSlices][i];
                sliceFieldRef[cellI] = nSlices_ + matSlices;
            }
            matSlices += 1;
        }

        // Check for empty slices
        forAll(sliceAddrList_, sliceI)
        {
            scalar sizeI(sliceAddrList_[sliceI].size());
            reduce(sizeI, sumOp<scalar>());

            if(sizeI <= 0)
            {
                FatalErrorIn("sliceMapperAutoAxialSlices")
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

Foam::sliceMapperAutoAxialSlices::sliceMapperAutoAxialSlices
(
    const fvMesh& mesh,
    const materials& mat,
    const dictionary& mapperOptDict
)
:
    sliceMapper(mesh, mat, mapperOptDict),
    precision_
    (       
#ifdef OPENFOAMFOUNDATION
        mapperOptDict.isNull() ? 
#elif OPENFOAMESI
        mapperOptDict.isNullDict() ? 
#endif             
        1e-6 :
        mapperOptDict.lookupOrDefault("precision", 1e-6)
    )
{   
    createSliceID();
    
    const globalOptions& globalOpt
    (mesh_.lookupObject<globalOptions>("globalOptions"));

    pinDirection_ = (globalOpt.pinDirection());

    calcAddressing();
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// ************************************************************************* //
