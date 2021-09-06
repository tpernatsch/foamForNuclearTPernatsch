/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

Description
    Heat exchanger.

Author
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

Class
    Foam::heatExchanger

SourceFiles
    heatExchanger.C

\*---------------------------------------------------------------------------*/

#ifndef heatExchanger_H
#define heatExchanger_H

#include "volFields.H"
#include "meshToMesh.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Forward declarations
class fluid;
class structure;

/*---------------------------------------------------------------------------*\
                           Class heatExchanger Declaration
\*---------------------------------------------------------------------------*/

class heatExchanger
: 
    public IOdictionary
{
protected :

    //- Ref to mesh
    const fvMesh& mesh_;

    //- Interfacial area, kinda useless to treat is as a volScalarField as
    //  unlike the structure and powerModels, there will be a separate
    //  heatExchanger object for each pair of cellZones
    scalar iA_;

    //- Wall conductance (i.e. wall thermal conductivity / wall thickness)
    scalar Hw_;

    //-
    labelList primaryCells_;

    //-
    labelList secondaryCells_;

    //-
    autoPtr<meshToMesh> mappingPtr_;

public:

    TypeName("heatExchanger");

    //- Constructors

        heatExchanger
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Destructor
        virtual ~heatExchanger() {}

    //- Member Functions

        //- Return volumetric area as a field rather than a scalar, so that it
        //  can also act as an HX marker field (i.e. iA = 0 in those cells where
        //  the HX does not exist)
        tmp<volScalarField> iA() const;

        //- In two-phase solvers :
        //  HT = H1*frac1*T1 + H2*frac2*T2
        //  H =  H1*frac1 + H2*frac2
        //- In a mono-phase solver :
        //  HT = H*T
        //  H = H
        //  THX is the field to be set according to the updated HX surface
        //  temperatures with might mapping magic
        void correct
        ( 
            const volScalarField& HT, 
            const volScalarField& H,
            volScalarField& THX
        );    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
