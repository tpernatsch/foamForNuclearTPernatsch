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
    Foam::PoissonRatioConstant

Description
    Model for constant Poisson Ratio. The value is read from dictionary.

SourceFiles
    PoissonRatioConstant.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef PoissonRatioConstant_H
#define PoissonRatioConstant_H

#include "PoissonRatioModel.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class PoissonRatioConstant Declaration
\*---------------------------------------------------------------------------*/

class PoissonRatioConstant
:
    public PoissonRatioModel
{
    // Private data
        
        //- Constant Poisson Ratio [-]
        dimensionedScalar nu_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        PoissonRatioConstant(const PoissonRatioConstant&);

        //- Disallow default bitwise assignment
        void operator=(const PoissonRatioConstant&);

protected:
    
    // Protected data

public:

    //- Runtime type information
        TypeName("constant");

    // Constructors

        //- Construct from mesh and dictionary
        PoissonRatioConstant
        (
            const fvMesh& mesh,
            const dictionary& dict, 
            const word defaultModel   
        );

    //- Destructor
        virtual ~PoissonRatioConstant();


    // Member Functions
    
    //- Update conductivity  
        virtual void correct
        (scalarField& sf, const scalarField& T, const labelList& addr);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
