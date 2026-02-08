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
    Foam::constantPoissonRatioUO2

Description
    Class to set UO2 Poisson's Ratio to a constant input value.

SourceFiles
    constantPoissonRatioUO2.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef constantPoissonRatioUO2_H
#define constantPoissonRatioUO2_H

#include "PoissonRatioModel.H"
#include "Switch.H"
#include "sliceMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class constantPoissonRatioUO2 Declaration
\*---------------------------------------------------------------------------*/

class constantPoissonRatioUO2
:
    public PoissonRatioModel
{
    // Private data

        //- Constant value 
        scalar PoissonRatioValue_;

        //- Activate isotopic cracking 
        Switch isotropicCracking_;

        //- Decide between alternative isotropic cracking models 
        word isotropicCrackingType_;

        //- Reference to nCracks field
        volScalarField& nCracks_;

        //- Maximum number of cracks (from 0 to 12)
        scalar nCracksMax_;
        
        //- Reference to slice mapper
        const sliceMapper* mapper_;   

    // Private Member Functions

        //- Disallow default bitwise copy construct
        constantPoissonRatioUO2(const constantPoissonRatioUO2&);

        //- Disallow default bitwise assignment
        void operator=(const constantPoissonRatioUO2&);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("UO2Constant");

    // Constructors

        //- Construct from mesh and dictionary
        constantPoissonRatioUO2
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dictionary& baseDict, 
            const word defaultModel   
        );

    //- Destructor
    virtual ~constantPoissonRatioUO2();


    // Member Functions
    
    //- Update conductivity  
    virtual void correct(scalarField& sf, const scalarField& T, const labelList& addr);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
