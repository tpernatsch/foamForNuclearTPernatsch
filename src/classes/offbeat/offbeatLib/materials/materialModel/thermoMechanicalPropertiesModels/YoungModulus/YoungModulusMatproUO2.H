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
    Foam::YoungModulusMatproUO2

Description
    Class modelling Young Modulus of UO2 fuel from MATPROv11.

SourceFiles
    YoungModulusMatproUO2.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef YoungModulusMatproUO2_H
#define YoungModulusMatproUO2_H

#include "YoungModulusModel.H"
#include "Switch.H"
#include "sliceMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class YoungModulusMatproUO2 Declaration
\*---------------------------------------------------------------------------*/

class YoungModulusMatproUO2
:
    public YoungModulusModel
{
    // Private data

        //- Fraction of Theoretical Density 
        scalar densityFrac_;

        //- Parameters for the Matpro model
        scalar par1;
        scalar par2;
        scalar par3;

        //- Perturbation parameter
        scalar perturb;

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
        YoungModulusMatproUO2(const YoungModulusMatproUO2&);

        //- Disallow default bitwise assignment
        void operator=(const YoungModulusMatproUO2&);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("UO2MATPRO");

    // Declare run-time constructor selection table

    // Constructors

        //- Construct from dictionary
        YoungModulusMatproUO2
        (
            const fvMesh& mesh,
            const dictionary& dict, 
            const word defaultModel   
        );

    //- Destructor
    virtual ~YoungModulusMatproUO2();


    // Member Functions
    
    //- Update conductivity  
    virtual void correct(scalarField& sf, const scalarField& T, const labelList& addr);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
