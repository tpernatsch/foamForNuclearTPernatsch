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
    Foam::thermalExpansionPARFUMESiC

Description
    Class modelling the thermal expansion of SiC with a constant value 4.9×10^−6 /K from PARFUME.

SourceFiles
    thermalExpansionPARFUMESiC.C

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

#ifndef thermalExpansionPARFUMESiC_H
#define thermalExpansionPARFUMESiC_H

#include "thermalExpansionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class thermalExpansionPARFUMESiC Declaration
\*---------------------------------------------------------------------------*/

class thermalExpansionPARFUMESiC
:
    public thermalExpansionModel
{
    // Private data

    //- The thermal expansion coefficient (1/K)
        scalar alpha_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        thermalExpansionPARFUMESiC(const thermalExpansionPARFUMESiC&);

        //- Disallow default bitwise assignment
        void operator=(const thermalExpansionPARFUMESiC&);

protected:

    // Protected data

public:

    //- Runtime type information
    TypeName("SiCPARFUME");

    // Constructors

        //- Construct from dictionary
        thermalExpansionPARFUMESiC
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word defaultModel
        );

    //- Destructor
    virtual ~thermalExpansionPARFUMESiC();


    // Member Functions

    //- Update thermalExpansion
    virtual void correct(symmTensorField& sf, const scalarField& T, const labelList& addr) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
