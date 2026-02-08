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
    Foam::sliceMapperByPellets

Description
    The slices are created with a similar algorithm as in the `byMaterial` mapper
    class. For fuel materials (e.g. UO2) the slice height coincide with the pellet
    size, introduced by the user in the material subdictionary. For non-fuel 
    materials, all cells are agglomerated into 1 slice.

    The algorithm tries to create one slice per pellet, which for 2D or 3D discrete 
    pellet simulations is probably the most logical approximation for phenomena such 
    as relocation or cracking (e.g. no variation of number of cracks along a pellet)

    For coarser simulations (e.g. in 1D or 2D smeared), the cells will be taller
    than the pellet size. The slice division algorithm will adapt following the 
    mesh axial discretization. For extruded meshes, this will create one slice
    per axial layer of cells.

    \warning
	This slice mapper class might duplicate functionalities that can be already
	achieved using the standard `byMaterial` class, so it might be eliminated
	in the future.

________________________________________________________________________________

\par Options
    The `byPellet` slice mapper requires only one option per fuel-type material in the 
    respective subdictionary located within the `solverDict`.

    Parameters in `materials` subdictionaries for fuel-type materials:

	- <b>`nPellets`</b> - Number of pellets in this material zone.

________________________________________________________________________________

Usage
    Here is a code snippet of the `solverDict` as an example of how to correctly 
    set the `byPellet` axial mapper:

    \verbatim
        ...
        sliceMapper     byPellet;

        materials
        {
            fuel
            {
                material   UO2;

                ...

                nPellets 100;

                ...
            }
        }
    \endverbatim 


SourceFiles
    sliceMapperByPellets.C

\todo
    -Check byMaterial slice mapper in multiprocessor
    -The addressing is based on the location of the cell center. Maybe other algorithms
    should be developed to allow for partial split of a cell between two consecutive
    slices? 
    -Possible introduction of cutting planes?

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)\n

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef sliceMapperByPellets_H
#define sliceMapperByPellets_H

#include "sliceMapper.H"
#include "globalOptions.H"
#include "fuelMaterial.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class sliceMapperByPellets Declaration
\*---------------------------------------------------------------------------*/

class sliceMapperByPellets
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
        sliceMapperByPellets(const sliceMapperByPellets&);

        //- Disallow default bitwise assignment
        void operator=(const sliceMapperByPellets&);

protected:
    
    // Protected data  


    // Protected member function

        // Update addressing sliceAddrList
        virtual void calcAddressing() const;
        
public:

    //- Runtime type information
    TypeName("byPellets");


    // Constructors

        //- Construct from mesh
        sliceMapperByPellets
        (
            const fvMesh& mesh,
            const materials& mat,
            const dictionary& sliceMapperDict
        );
    

    //- Destructor
    ~sliceMapperByPellets()
    {}


    // Member Functions    

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
