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
    Foam::YoungModulusSckCenUPuO2

Description
    Correlation for Young's modulus of (U,Pu)O2 MOX fuel from SCK-CEN
    developed for TRANSURANUS.

    Source: https://re.public.polimi.it/retrieve/handle/11311/1172415/619939/Technical%20Report_INSPYRE_WP7-D7.2_%282020%29.pdf

SourceFiles
    YoungModulusSckCenUPuO2.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef YoungModulusSckCenUPuO2_H
#define YoungModulusSckCenUPuO2_H

#include "YoungModulusModel.H"
#include "Switch.H"
#include "sliceMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class YoungModulusSckCenUPuO2 Declaration
\*---------------------------------------------------------------------------*/

class YoungModulusSckCenUPuO2
:
    public YoungModulusModel
{
    // Private data

        //- Deviation from stoichiometry y = 2-O/M
        scalar y_;

        //- Pu/HM fraction
        scalar ratioPuMetal_;

        //- Model Parameters
        scalar dEdy_;
        scalar E1_;
        scalar E2_;
        scalar E3_;

        //- Young's moduli UO2 UPuO2
        scalar EUO2_;
        scalar EPuO2_;

        //- Perturbation parameter
        scalar perturb;

        //- Activate isotopic cracking 
        Switch isotropicCracking_;

        //- Reference to nCracks field
        volScalarField& nCracks_;

        //- Maximum number of cracks (from 0 to 12)
        scalar nCracksMax_;
    
        //- Reference to slice mapper
        const sliceMapper* mapper_;   

    // Private Member Functions

        //- Disallow default bitwise copy construct
        YoungModulusSckCenUPuO2(const YoungModulusSckCenUPuO2&);

        //- Disallow default bitwise assignment
        void operator=(const YoungModulusSckCenUPuO2&);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("UPuO2SckCen");

    // Declare run-time constructor selection table

    // Constructors

        //- Construct from dictionary
        YoungModulusSckCenUPuO2
        (
            const fvMesh& mesh,
            const dictionary& dict, 
            const word defaultModel   
        );

    //- Destructor
    virtual ~YoungModulusSckCenUPuO2();


    // Member Functions
    
        //- Update conductivity  
        virtual void correct(scalarField& sf, const scalarField& T, const labelList& addr);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
