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
    Foam::ZircaloyOuterCorrosion

Description
    Water-side corrosion for Zircaloy cladding.
	Currently hardcoded for the EPRI/KWU/C-E low-temperature oxidation
	model and Cathcart-Pawel high temperature oxidation model.
	Model constants are taken from (Dunbar et al, 2024).

	Dunbar et al, Fuel performance analysis of Cr-coated Zircaloy-4 cladding
	during a prototypical LOCA event using BISON, Annals of Nuclear Energy
	Volume 200, 1 June 2024, https://doi.org/10.1016/j.anucene.2024.110411.

Usage
    See Foam::corrosionByPatch

SourceFiles
    ZircaloyOuterCorrosion.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef ZircaloyOuterCorrosion_H
#define ZircaloyOuterCorrosion_H

#include "corrosionModel.H"
#include "oxidationKineticsModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class ZircaloyOuterCorrosion Declaration
\*---------------------------------------------------------------------------*/

class ZircaloyOuterCorrosion
:
    public corrosionModel
{
    // Private data

        autoPtr<oxidationKineticsModel> kineticsModel_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        ZircaloyOuterCorrosion(const ZircaloyOuterCorrosion&);

        //- Disallow default bitwise assignment
        void operator=(const ZircaloyOuterCorrosion&);

protected:

    // Protected data

    // Protected Member Functions

public:

    //- Runtime type information
    TypeName("ZircaloyOuterCorrosion");


    // Constructors

        //- Construct from dictionary
        ZircaloyOuterCorrosion
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Destructors

        virtual ~ZircaloyOuterCorrosion();


    // Access functions


    // Member Functions

        //- Update oxide thickness
        virtual void correct
        (
            surfaceScalarField& oxideThickness,
            surfaceScalarField& DOxideThickness,
            const label& patchID
        );

        //- Update change in metal thickness
        virtual void updateDMetalThickness
        (
            const surfaceScalarField& DOxideThickness,
            surfaceScalarField& DMetalThickness,
            const label& patchID
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
