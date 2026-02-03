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
    Foam::conductivityUPuO2LanningBeyer

Description

SourceFiles
    conductivityUPuO2LanningBeyer.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef conductivityUPuO2LanningBeyer_H
#define conductivityUPuO2LanningBeyer_H

#include "conductivityModel.H"
#include "materialModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class conductivityUPuO2LanningBeyer Declaration
\*---------------------------------------------------------------------------*/

class conductivityUPuO2LanningBeyer
:
    public conductivityModel
{
    // Private data

        //- Name used for the burnup field
        const word burnupName_;
        
        //- Reference to Burnup field
        const volScalarField* Bu_;

        //- Thermal conductivity for a 95% TD MOX fuel [W/m K]
        scalar kappa95_;

        //- O/M - oxygen to metal ratio
        const volScalarField& OM_;

        // Fuel density fraction
        scalar densityFrac_;

        // Fuel porosity
        scalar porosity_;
        
    // Private Member Functions

        //- Disallow default bitwise copy construct
        conductivityUPuO2LanningBeyer(const conductivityUPuO2LanningBeyer&);

        //- Disallow default bitwise assignment
        void operator=(const conductivityUPuO2LanningBeyer&);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("UPuO2LanningBeyer");

    // Constructors

        //- Construct from dictionary
        conductivityUPuO2LanningBeyer
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dictionary& baseDict,
            const word defaultModel         
        );
    
    //- Destructor
    virtual ~conductivityUPuO2LanningBeyer();


    // Member Functions
    
        //- Update conductivity  
        virtual void correct(scalarField& sf, const scalarField& T, const labelList& addr);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
