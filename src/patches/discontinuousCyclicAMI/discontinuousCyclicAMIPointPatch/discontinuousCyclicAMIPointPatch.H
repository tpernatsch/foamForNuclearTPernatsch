/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::cyclicAMIPointPatch

Description
    Cyclic AMI point patch - place holder only

SourceFiles
    cyclicAMIPointPatch.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_discontinuousCyclicAMIPointPatch_H
#define Foam_discontinuousCyclicAMIPointPatch_H

#include "cyclicAMIPointPatch.H"

namespace Foam
{

class discontinuousCyclicAMIPointPatch
:
    public cyclicAMIPointPatch
{
    // Private Member Functions

        //- No copy construct
        discontinuousCyclicAMIPointPatch(const discontinuousCyclicAMIPointPatch&) = delete;

        //- No copy assignment
        void operator=(const discontinuousCyclicAMIPointPatch&) = delete;

public:

    //- Runtime type information
    TypeName("discontinuousCyclicAMI");


    // Constructors

        //- Construct from components
        discontinuousCyclicAMIPointPatch
        (
            const polyPatch& patch,
            const pointBoundaryMesh& bm
        );


    //- Destructor
    virtual ~discontinuousCyclicAMIPointPatch();


    // Member Functions

        //- Return the constraint type this pointPatch implements.
        //  For discontinuous, this is disabled.
        virtual const word& constraintType() const
        {
            return word::null;
        }
};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //