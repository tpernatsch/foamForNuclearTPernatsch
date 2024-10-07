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
    Foam::YoungModulusPARFUMESiC

Description
    Class modelling Young Modulus of SiC with the Interpolation from PARFUME.

SourceFiles
    YoungModulusPARFUMESiC.C

    \htmlinclude YoungModulus_YoungModulusPARFUMESiC.html

\mainauthor
    F. Xiang - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
    Laboratory for Reactor Physics and Systems Behaviour) & XJTU (Xi'an JiaoTong
    University, China, Nuclear Thermal-hydraulic Laboratory)

\contribution
    A. Scolaro - EPFL
    E. Brunetto - EPFL

\date
    December 2022

\*---------------------------------------------------------------------------*/

#ifndef YoungModulusPARFUMESiC_H
#define YoungModulusPARFUMESiC_H

#include "YoungModulusModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class YoungModulusPARFUMESiC Declaration
\*---------------------------------------------------------------------------*/

class YoungModulusPARFUMESiC
:
    public YoungModulusModel
{
    // Private data;

        //- Parameters for the Interpolation
        scalar temp1;
        scalar temp2;
        scalar temp3;
        scalar temp4;
        scalar E1;
        scalar E2;
        scalar E3;
        scalar E4;

        //- Perturbation parameter
        scalar perturb;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        YoungModulusPARFUMESiC(const YoungModulusPARFUMESiC&);

        //- Disallow default bitwise assignment
        void operator=(const YoungModulusPARFUMESiC&);

protected:

    // Protected data

public:

    //- Runtime type information
    TypeName("SiCPARFUME");

    // Constructors

        //- Construct from dictionary
        YoungModulusPARFUMESiC
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word defaultModel
        );
    //- Destructor
    virtual ~YoungModulusPARFUMESiC();


    // Member Functions

    //- Update conductivity
    virtual void correct(scalarField& sf, const scalarField& T, const labelList& addr) ;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
