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
    Foam::heatCapacitySneadSiC

Description
    Correlation for SiC heat capacity from Snead et al. (2007).

    \htmlinclude heatCapacity_heatCapacitySneadSiC.html

SourceFiles
    heatCapacitySneadSiC.C

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

#ifndef heatCapacitySneadSiC_H
#define heatCapacitySneadSiC_H

#include "heatCapacityModel.H"
#include "physicoChemicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class heatCapacitySneadSiC Declaration
\*---------------------------------------------------------------------------*/

class heatCapacitySneadSiC
:
    public heatCapacityModel
{
    // Private data

        //- Correlation parameters
        scalar par1;
        scalar par2;
        scalar par3;
        scalar par4;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        heatCapacitySneadSiC(const heatCapacitySneadSiC&);

        //- Disallow default bitwise assignment
        void operator=(const heatCapacitySneadSiC&);

protected:

    // Protected data

public:

    //- Runtime type information
    TypeName("SiCSnead");

    // Declare run-time constructor selection table

    // Constructors

        //- Construct from mesh, dicionary and defaultModel
        heatCapacitySneadSiC
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dictionary& baseDict,
            const word defaultModel
        );


    // Selectors


    //- Destructor
    virtual ~heatCapacitySneadSiC();


    // Member Functions

    //- Update heatCapacity
    virtual void correct(scalarField& sf, const scalarField& T, const labelList& addr) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
