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
    Foam::constantPoissonRatioMolybdenum

Description
    Class to set Mo Poisson's Ratio to a constant input value.

SourceFiles
    constantPoissonRatioMolybdenum.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef constantPoissonRatioMolybdenumUO2_H
#define constantPoissonRatioMolybdenum_H

#include "PoissonRatioModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class constantPoissonRatioMolybdenum Declaration
\*---------------------------------------------------------------------------*/

class constantPoissonRatioMolybdenum
:
    public PoissonRatioModel
{
    // Private data

        //- Constant value 
        scalar PoissonRatioValue_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        constantPoissonRatioMolybdenum(const constantPoissonRatioMolybdenum&);

        //- Disallow default bitwise assignment
        void operator=(const constantPoissonRatioMolybdenum&);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("MoConstant");

    // Constructors

        //- Construct from mesh and dictionary
        constantPoissonRatioMolybdenum
        (
            const fvMesh& mesh,
            const dictionary& dict, 
            const word defaultModel   
        );

    //- Destructor
    virtual ~constantPoissonRatioMolybdenum();


    // Member Functions
    
    //- Update conductivity  
    virtual void correct(scalarField& sf, const scalarField& T, const labelList& addr);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
