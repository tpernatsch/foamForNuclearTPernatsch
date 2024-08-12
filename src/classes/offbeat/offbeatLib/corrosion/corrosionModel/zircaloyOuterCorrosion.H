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
    Foam::zircaloyOuterCorrosion

Description
    Water-side corrosion for zircaloy cladding. The model's constants are taken
    from the BISON manual.

Usage

SourceFiles
    zircaloyOuterCorrosion.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef zircaloyOuterCorrosion_H
#define zircaloyOuterCorrosion_H

#include "corrosionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class zircaloyOuterCorrosion Declaration
\*---------------------------------------------------------------------------*/

class zircaloyOuterCorrosion
:
    public corrosionModel
{
    // Private data
    
    // Private Member Functions

        //- Disallow default bitwise copy construct
        zircaloyOuterCorrosion(const zircaloyOuterCorrosion&);

        //- Disallow default bitwise assignment
        void operator=(const zircaloyOuterCorrosion&);

protected:
    
    // Protected data
    
    // Protected Member Functions

public:

    //- Runtime type information
    TypeName("zircaloyOuterCorrosion");


    // Constructors

        //- Construct from dictionary
        zircaloyOuterCorrosion
        (
            const fvMesh& mesh
        );


    // Selectors


    // Destructors

        virtual ~zircaloyOuterCorrosion();


    // Access functions
        

    // Member Functions
    
        //- Update oxide thickness 
        virtual void correct
        (
            surfaceScalarField& oxideThickness,
            surfaceScalarField& DOxideThickness,
            const label& patchID
        );
    
        //- Update change in metal thickness 
        virtual void updateDMetalThickness
        (
            const surfaceScalarField& DOxideThickness,
            surfaceScalarField& DMetalThickness,
            const label& patchID
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
