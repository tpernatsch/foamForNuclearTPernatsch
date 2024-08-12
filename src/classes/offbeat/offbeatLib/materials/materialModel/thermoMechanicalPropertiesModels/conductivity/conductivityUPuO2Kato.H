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
    Foam::conductivityUPuO2Kato

Description
    S.Novascone et al.,Modeling porosity migration in LWR and fast reactor MOX 
    fuel using the finite element method, Journal of Nuclear Materials, Volume 
    508, 2018, Pages 226-236, ISSN 0022-3115, 
    https://doi.org/10.1016/j.jnucmat.2018.05.041.
    (https://www.sciencedirect.com/science/article/pii/S0022311518302861)

SourceFiles
    conductivityUPuO2Kato.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef conductivityUPuO2Kato_H
#define conductivityUPuO2Kato_H

#include "conductivityModel.H"
#include "materialModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class conductivityUPuO2Kato Declaration
\*---------------------------------------------------------------------------*/

class conductivityUPuO2Kato
:
    public conductivityModel
{
    // Private data

        //- Initial density fraction
        const scalar densityFrac_;

        //- Ptr to porosity internal field
        const volScalarField* p_;

        //- Am weight fraction
        const scalar Am_;

        //- Np weight fraction
        const scalar Np_;
        
        //- Stoichiometry deviation
        const scalar x_;

        //- Perturbation factor
        scalar perturb;
    
    // Private Member Functions

        //- Disallow default bitwise copy construct
        conductivityUPuO2Kato(const conductivityUPuO2Kato&);

        //- Disallow default bitwise assignment
        void operator=(const conductivityUPuO2Kato&);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("UPuO2Kato");

    // Constructors

        //- Construct from dictionary
        conductivityUPuO2Kato
        (
            const fvMesh& mesh,
            const dictionary& dict, 
            const word defaultModel   
        );
    
    //- Destructor
    virtual ~conductivityUPuO2Kato();


    // Member Functions
    
        //- Update conductivity  
        virtual void correct(scalarField& sf, const scalarField& T, const labelList& addr);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
