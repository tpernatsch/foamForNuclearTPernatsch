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
    Foam::conductivityUPuO2Brancheria

Description
    Conductivity correlation for MOX fuel Pu=0.25, U=0.75 from Brancheria et al.
    Ref : "Calibration of a fuel-to-cladding gap conductance model for fast 
    reactor fuel pins", Hanford Engineering Development Laboratory 
    url : https://www.osti.gov/servlets/purl/6863540

SourceFiles
    conductivityUPuO2Brancheria.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro, C. Fiorina - EPFL\n

\date 
    February 2024

\*---------------------------------------------------------------------------*/

#ifndef conductivityUPuO2Brancheria_H
#define conductivityUPuO2Brancheria_H

#include "conductivityModel.H"
#include "materialModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class conductivityUPuO2Brancheria Declaration
\*---------------------------------------------------------------------------*/

class conductivityUPuO2Brancheria
:
    public conductivityModel
{
    // Private data

        //- Initial density fraction
        const scalar densityFrac_;

        //- Ptr to porosity internal field
        const volScalarField* p_;

        //- Perturbation factor
        scalar perturb;
    
    // Private Member Functions

        //- Disallow default bitwise copy construct
        conductivityUPuO2Brancheria(const conductivityUPuO2Brancheria&);

        //- Disallow default bitwise assignment
        void operator=(const conductivityUPuO2Brancheria&);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("UPuO2Brancheria");

    // Constructors

        //- Construct from dictionary
        conductivityUPuO2Brancheria
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dictionary& baseDict,
            const word defaultModel         
        );
    
    //- Destructor
    virtual ~conductivityUPuO2Brancheria();


    // Member Functions
    
        //- Update conductivity  
        virtual void correct(scalarField& sf, const scalarField& T, const labelList& addr);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
