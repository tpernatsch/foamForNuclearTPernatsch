/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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
    Foam::constantEmissivityMolybdenum

Description
    Class to set Mo emissivity to a constant input value.

SourceFiles
    constantEmissivityMolybdenum.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef constantEmissivityMolybdenum_H
#define constantEmissivityMolybdenum_H

#include "emissivityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class constantEmissivityMolybdenum Declaration
\*---------------------------------------------------------------------------*/

class constantEmissivityMolybdenum
:
    public emissivityModel
{
    // Private data

        //- Constant emissitvity value
        scalar emissivityValue_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        constantEmissivityMolybdenum(const constantEmissivityMolybdenum&);

        //- Disallow default bitwise assignment
        void operator=(const constantEmissivityMolybdenum&);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("MoConstant");

    // Constructors

        //- Construct from dictionary
        constantEmissivityMolybdenum
        (
            const fvMesh& mesh,
            const dictionary& dict, 
            const word defaultModel   
        );

    //- Destructor
    virtual ~constantEmissivityMolybdenum();


    // Member Functions
    
    //- Update emissivity  
    virtual void correct(scalarField& sf, const scalarField& T, const labelList& addr) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
