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
    Foam::sliceMapperByMaterial

Description
    The slices are created on a material basis. The number and height of the 
    slices are provided by the user in the materials subdictionary.

Usage
    In solverDict file:
    \verbatim
    ...
    sliceMapper     byMaterial;

    materials
    {
        fuel
        {
            material   UO2;

            ...

            nSlices 10;

            // If heightSlices is not specified, the fuel will be split
            // nSlices uniform slices.
            // In this example, the fuel (0.5m long) is split in 5 slices of 
            // 0.1m
            heightSlices( 0.1 0.1 0.1 0.1 0.1);

            ...
        }
    }
    \endverbatim 

\todo
    -Check byMaterial slice mapper in multiprocessor
    -The addressing is based on the location of the cell center. Maybe other algorithms
    should be developed to allow for partial split of a cell between two consecutive
    slices? 
    -Possible introduction of cutting planes?

SourceFiles
    sliceMapperByMaterial.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)\n

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef sliceMapperByMaterial_H
#define sliceMapperByMaterial_H

#include "sliceMapper.H"
#include "globalOptions.H"
#include "fuelMaterial.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class sliceMapperByMaterial Declaration
\*---------------------------------------------------------------------------*/

class sliceMapperByMaterial
:
    public sliceMapper
{
    // Private data

        //- Direction of fuel pin central axis
        vector pinDirection_;
    

    // Private Member Functions
        
        //- Check that all cells are accounted for in zones
        void checkZoneDefinition();

        //- Disallow default bitwise copy construct
        sliceMapperByMaterial(const sliceMapperByMaterial&);

        //- Disallow default bitwise assignment
        void operator=(const sliceMapperByMaterial&);

protected:
    
    // Protected data  


    // Protected member function

        // Update addressing sliceAddrList
        virtual void calcAddressing() const;
        
public:

    //- Runtime type information
    TypeName("byMaterial");


    // Constructors

        //- Construct from mesh
        sliceMapperByMaterial
        (
            const fvMesh& mesh,
            const materials& mat,
            const dictionary& mapperOptDict
        );
    

    //- Destructor
    ~sliceMapperByMaterial()
    {}


    // Member Functions
    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
