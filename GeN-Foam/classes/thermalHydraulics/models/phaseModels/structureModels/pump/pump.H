/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

Description
    Pump.

Author
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

Class
    Foam::pump

SourceFiles
    pump.C

\*---------------------------------------------------------------------------*/

#ifndef pump_H
#define pump_H

#include "volFields.H"
#include "meshToMesh.H"
#include "zeroGradientFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Forward declarations
class structureModel;

/*---------------------------------------------------------------------------*\
                           Class pump Declaration
\*---------------------------------------------------------------------------*/

class pump
: 
    public IOdictionary
{
protected :

    //- Ref to mesh
    const fvMesh& mesh_;

    //- List of cells in which this pump exists
    const labelList& cellList_;

    //- Const value of the momentum source
    const vector pumpValue_;
    
    //- Optional table to describe any possible dependence of the momentum
    //  source on time
    autoPtr<Function1<scalar>> timeProfilePtr_;

    //- Optional offset time
    scalar t0_;

    //- True if a timeProfilePtr_ exists
    bool timeDependent_;

public:

    TypeName("pump");

    //- Constructors

        pump
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const labelList& cellList
        );

        //- Destructor
        virtual ~pump() {}

    //- Member Functions

        //- Update momentum source field
        void correct(volVectorField& pumpField);    

        //- Access

            //-
            bool timeDependent() const
            {
                return timeDependent_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
