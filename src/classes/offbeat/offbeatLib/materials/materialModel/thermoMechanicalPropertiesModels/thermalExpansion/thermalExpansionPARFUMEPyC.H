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
    Foam::thermalExpansionPARFUMEPyC

Description
    Class modelling the transversly isotropic thermal expansion of PyC through
    the code PARFUME.

    \htmlinclude thermalExpansion_thermalExpansionPARFUMEPyC.html

SourceFiles
    thermalExpansionPARFUMEPyC.C

\mainauthor
    F. Xiang - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
    Laboratory for Reactor Physics and Systems Behaviour) & XJTU (Xi'an JiaoTong
    University, China, Nuclear Thermal-hydraulic Laboratory)

\contribution
    A. Scolaro - EPFL
    E. Brunetto - EPFL

\date
    November 2022

\*---------------------------------------------------------------------------*/

#ifndef thermalExpansionPARFUMEPyC_H
#define thermalExpansionPARFUMEPyC_H

#include "thermalExpansionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class thermalExpansionPARFUMEPyC Declaration
\*---------------------------------------------------------------------------*/

class thermalExpansionPARFUMEPyC
:
    public thermalExpansionModel
{
    // Private data

        //- The as-fabricated anisotropy (dimensionless)
        scalar BAF_;

        //- Whether to use spherical coordinate
        bool sp_;

        //- Parameters for the PARFUME model
        scalar par1;
        scalar par2;
        scalar par3;
        scalar par4;
        scalar par5;
        scalar par6;

        //- Perturbation parameter
        scalar perturb;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        thermalExpansionPARFUMEPyC(const thermalExpansionPARFUMEPyC&);

        //- Disallow default bitwise assignment
        void operator=(const thermalExpansionPARFUMEPyC&);

protected:

    // Protected data

public:

    //- Runtime type information
    TypeName("PyCPARFUME");

    // Constructors

        //- Construct from dictionary
        thermalExpansionPARFUMEPyC
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word defaultModel
        );

    //- Destructor
    virtual ~thermalExpansionPARFUMEPyC();


    // Member Functions

    //- Update thermalExpansion
    virtual void correct(symmTensorField& sf, const scalarField& T, const labelList& addr) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
