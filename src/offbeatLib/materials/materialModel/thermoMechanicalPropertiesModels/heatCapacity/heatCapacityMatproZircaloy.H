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
    Foam::heatCapacityMatproZircaloy

Description
    Correlation for Zircaloy heat capacity from Matpro.
    Link : https://www.nrc.gov/docs/ML1429/ML14296A063.pdf - page 60.

SourceFiles
    heatCapacityMatproZircaloy.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef heatCapacityMatproZircaloy_H
#define heatCapacityMatproZircaloy_H

#include "heatCapacityModel.H"
#include "physicoChemicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class heatCapacityMatproZircaloy Declaration
\*---------------------------------------------------------------------------*/

class heatCapacityMatproZircaloy
:
    public heatCapacityModel
{
    // Private data
    
    // Private Member Functions

        //- Disallow default bitwise copy construct
        heatCapacityMatproZircaloy(const heatCapacityMatproZircaloy&);

        //- Disallow default bitwise assignment
        void operator=(const heatCapacityMatproZircaloy&);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("ZircaloyMatpro");

    // Declare run-time constructor selection table

    // Constructors

        //- Construct from mesh, dicionary and defaultModel
        heatCapacityMatproZircaloy
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dictionary& baseDict, 
            const word defaultModel
        );


    // Selectors


    //- Destructor
    virtual ~heatCapacityMatproZircaloy();


    // Member Functions
    
    //- Update heatCapacity  
    virtual void correct(scalarField& sf, const scalarField& T, const labelList& addr) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
