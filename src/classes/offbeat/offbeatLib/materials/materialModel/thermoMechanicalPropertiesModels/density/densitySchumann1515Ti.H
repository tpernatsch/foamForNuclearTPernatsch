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
    Foam::densitySchumann1515Ti

Description
    Class modelling the density evolution of 15-15Ti cladding material through 
    the Schumann (1970) correlation.

SourceFiles
    densitySchumann1515Ti.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef densitySchumann1515Ti_H
#define densitySchumann1515Ti_H

#include "densityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class densitySchumann1515Ti Declaration
\*---------------------------------------------------------------------------*/

class densitySchumann1515Ti
:
    public densityModel
{
    // Private data

        //- Model parameters
        scalar rho0_;

        //- Perturbation factor
        scalar perturb;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        densitySchumann1515Ti(const densitySchumann1515Ti&);

        //- Disallow default bitwise assignment
        void operator=(const densitySchumann1515Ti&);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("1515TiSchumann");

    // Constructors

        //- Construct from dictionary
        densitySchumann1515Ti
        (
            const fvMesh& mesh,
            const dictionary& dict, 
            const word defaultModel   
        );

    //- Destructor
    virtual ~densitySchumann1515Ti();


    // Member Functions
    
    //- Update conductivity  
    virtual void correct(scalarField& sf, const scalarField& T, const labelList& addr) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
