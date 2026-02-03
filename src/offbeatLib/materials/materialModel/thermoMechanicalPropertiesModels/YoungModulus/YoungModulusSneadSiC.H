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
    Foam::YoungModulusSneadSiC

Description
    Class modelling Young Modulus of SiC from Snead et al. work:
    L. L. Snead, T. Nozawa, Y. Katoh, T.-S. Byun, S. Kondo, and D. A. Petti.
    Handbook of sic properties for fuel performance modeling.
    Journal of Nuclear Materials, 371:329–377, 2007.

SourceFiles
    YoungModulusSneadSiC.C

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

#ifndef YoungModulusSneadSiC_H
#define YoungModulusSneadSiC_H

#include "YoungModulusModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class YoungModulusSneadSiC Declaration
\*---------------------------------------------------------------------------*/

class YoungModulusSneadSiC
:
    public YoungModulusModel
{
    // Private data;

        //- SiC porosity
        scalar Vp_;

        //- Parameters for the Snead model
        scalar E0;
        scalar T0;
        scalar B;
        scalar C;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        YoungModulusSneadSiC(const YoungModulusSneadSiC&);

        //- Disallow default bitwise assignment
        void operator=(const YoungModulusSneadSiC&);

protected:

    // Protected data

public:

    //- Runtime type information
    TypeName("SiCSnead");

    // Constructors

        //- Construct from dictionary
        YoungModulusSneadSiC
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dictionary& baseDict,
            const word defaultModel
        );
    //- Destructor
    virtual ~YoungModulusSneadSiC();


    // Member Functions

    //- Update conductivity
    virtual void correct(scalarField& sf, const scalarField& T, const labelList& addr) ;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
