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
    Foam::constantPoissonRatioUPuO2

Description
    Class modelling the constant Poisson's ratio of UPuO2 MOX fuel from 
    Matprov11.

SourceFiles
    constantPoissonRatioUPuO2.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef constantPoissonRatioUPuO2_H
#define constantPoissonRatioUPuO2_H

#include "PoissonRatioModel.H"
#include "sliceMapper.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class constantPoissonRatioUPuO2 Declaration
\*---------------------------------------------------------------------------*/

class constantPoissonRatioUPuO2
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
        constantPoissonRatioUPuO2(const constantPoissonRatioUPuO2&);

        //- Disallow default bitwise assignment
        void operator=(const constantPoissonRatioUPuO2&);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("UPuO2Constant");

    // Constructors

        //- Construct from mesh and dictionary
        constantPoissonRatioUPuO2
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dictionary& baseDict, 
            const word defaultModel   
        );

    //- Destructor
    virtual ~constantPoissonRatioUPuO2();


    // Member Functions
    
    //- Update conductivity  
    virtual void correct(scalarField& sf, const scalarField& T, const labelList& addr);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
