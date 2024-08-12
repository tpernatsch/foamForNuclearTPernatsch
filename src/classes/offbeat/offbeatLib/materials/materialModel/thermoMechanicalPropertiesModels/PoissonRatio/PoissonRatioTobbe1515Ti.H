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
    Foam::PoissonRatioTobbe1515Ti

Description
    Class modelling Poisson's ratio of 15-15 Ti cladding material based on Tobbe
    correlation (1975).

SourceFiles
    PoissonRatioTobbe1515Ti.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef PoissonRatioTobbe1515Ti_H
#define PoissonRatioTobbe1515Ti_H

#include "PoissonRatioModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class PoissonRatioTobbe1515Ti Declaration
\*---------------------------------------------------------------------------*/

class PoissonRatioTobbe1515Ti
:
    public PoissonRatioModel
{
    // Private data

        //- Parameters for the Matpro model
        scalar par1;
        scalar par2;

        //- Perturbation parameter
        scalar perturb;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        PoissonRatioTobbe1515Ti(const PoissonRatioTobbe1515Ti&);

        //- Disallow default bitwise assignment
        void operator=(const PoissonRatioTobbe1515Ti&);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("1515TiTobbe");

    // Constructors

        //- Construct from mesh and dictionary
        PoissonRatioTobbe1515Ti
        (
            const fvMesh& mesh,
            const dictionary& dict, 
            const word defaultModel   
        );

    //- Destructor
    virtual ~PoissonRatioTobbe1515Ti();


    // Member Functions
    
    //- Update conductivity  
    virtual void correct(scalarField& sf, const scalarField& T, const labelList& addr);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
