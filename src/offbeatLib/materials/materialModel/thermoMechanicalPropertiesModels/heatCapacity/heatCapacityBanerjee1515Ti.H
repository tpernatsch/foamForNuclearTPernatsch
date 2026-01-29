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
    Foam::heatCapacityBanerjee1515Ti

Description
    Correlation for 15-15 Ti heat capacity from Banerjee et al. (2007).

SourceFiles
    heatCapacityBanerjee1515Ti.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef heatCapacityBanerjee1515Ti_H
#define heatCapacityBanerjee1515Ti_H

#include "heatCapacityModel.H"
#include "physicoChemicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class heatCapacityBanerjee1515Ti Declaration
\*---------------------------------------------------------------------------*/

class heatCapacityBanerjee1515Ti
:
    public heatCapacityModel
{
    // Private data

        //- Correlation parameters
        scalar par1;
        scalar par2;
        scalar par3;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        heatCapacityBanerjee1515Ti(const heatCapacityBanerjee1515Ti&);

        //- Disallow default bitwise assignment
        void operator=(const heatCapacityBanerjee1515Ti&);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("Steel1515TiBanerjee");

    // Declare run-time constructor selection table

    // Constructors

        //- Construct from mesh, dicionary and defaultModel
        heatCapacityBanerjee1515Ti
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dictionary& baseDict, 
            const word defaultModel
        );


    // Selectors


    //- Destructor
    virtual ~heatCapacityBanerjee1515Ti();


    // Member Functions
    
    //- Update heatCapacity  
    virtual void correct(scalarField& sf, const scalarField& T, const labelList& addr) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
