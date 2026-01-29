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
    Foam::thermalExpansionMAMOX

Description
    Class modelling  isotropic thermal expansion in MA-MOX fuel (for Pu = 0.3) 
    Source : JNM 469 (2016) 223-227.
    url : https://doi.org/10.1016/j.jnucmat.2015.11.048

SourceFiles
    thermalExpansionMAMOX.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef thermalExpansionMAMOX_H
#define thermalExpansionMAMOX_H

#include "thermalExpansionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class thermalExpansionMAMOX Declaration
\*---------------------------------------------------------------------------*/

class thermalExpansionMAMOX
:
    public thermalExpansionModel
{
    // Private data

        //- Pu conc
        scalar cPu_;

        //- Oxygen to metal ratio
        scalar OM_;

        //- Deviation from stoichiometry x
        scalar x_;

        //- Perturbation parameter
        scalar perturb;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        thermalExpansionMAMOX(const thermalExpansionMAMOX&);

        //- Disallow default bitwise assignment
        void operator=(const thermalExpansionMAMOX&);

        //- Function to compute a0
        virtual const scalar a0( scalar cPu, scalar x) const;
         
        //- Function to compute a1
        virtual const scalar a1( scalar cPu, scalar x) const;
         
        //- Function to compute a2
        virtual const scalar a2( scalar cPu, scalar x) const;
         
        //- Function to compute a3
        virtual const scalar a3( scalar cPu, scalar x) const;
         


protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("MaUPuO2Kato");

    // Constructors

        //- Construct from dictionary
        thermalExpansionMAMOX
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dictionary& baseDict, 
            const word defaultModel   
        );

    //- Destructor
    virtual ~thermalExpansionMAMOX();


    // Member Functions
    
        //- Update thermalExpansion  
        virtual void correct(symmTensorField& sf, const scalarField& T, const labelList& addr) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
