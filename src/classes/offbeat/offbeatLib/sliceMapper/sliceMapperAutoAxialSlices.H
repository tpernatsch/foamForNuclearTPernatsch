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

Class
    Foam::sliceMapperAutoAxialSlices

Description
    Slice mapper developed for quasi2D or 2D r-z coarse smeared geometries.

    The slices are automaticaly created for each material.
    First, the cell centers and IDs are saved in a HashTable according to their
    axial location (i.e. the HashTable will have as many entries as there are 
    different axial locations in the material mesh). 
    Then the ID lists for each axial location in the HashTable are passed in 
    the sliceAddrList pointer list, making sure that the labelLists are passed in
    ascending axial location order (from bottom to top).


Usage
    In solverDict file:
    \verbatim

    sliceMapper     autoAxialSlices;
    
    \endverbatim 

\todo
    -Check byMaterial slice mapper in multiprocessor
    -The addressing is based on the location of the cell center. Maybe other algorithms
    should be developed to allow for partial split of a cell between two consecutive
    slices? 
    -Possible introduction of cutting planes?

SourceFiles
    sliceMapperAutoAxialSlices.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)\n

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef sliceMapperByPellet_H
#define sliceMapperByPellet_H

#include "sliceMapper.H"
#include "globalOptions.H"
#include "fuelMaterial.H"
#include <set>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class sliceMapperAutoAxialSlices Declaration
\*---------------------------------------------------------------------------*/

class sliceMapperAutoAxialSlices
:
    public sliceMapper
{
    // Private data

        //- Direction of fuel pin central axis
        vector pinDirection_;

        //- Precision of automatic slice division algorithm
        scalar precision_;
    

    // Private Member Functions
        
        //- Check that all cells are accounted for in zones
        void checkZoneDefinition();

        //- Disallow default bitwise copy construct
        sliceMapperAutoAxialSlices(const sliceMapperAutoAxialSlices&);

        //- Disallow default bitwise assignment
        void operator=(const sliceMapperAutoAxialSlices&);

protected:
    
    // Protected data  


    // Protected member function

        // Update addressing sliceAddrList
        virtual void calcAddressing() const;
        
public:

    //- Runtime type information
    TypeName("autoAxialSlices");


    // Constructors

        //- Construct from mesh
        sliceMapperAutoAxialSlices
        (
            const fvMesh& mesh,
            const materials& mat,
            const dictionary& mapperOptDict
        );
    

    //- Destructor
    ~sliceMapperAutoAxialSlices()
    {}


    // Member Functions    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
