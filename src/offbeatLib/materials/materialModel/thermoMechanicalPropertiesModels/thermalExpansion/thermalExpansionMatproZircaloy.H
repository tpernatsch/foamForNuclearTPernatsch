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
    Foam::thermalExpansionMatproZircaloy

Description
    Class modelling thermal expansion of UO2 fuel from Matprov11.

SourceFiles
    thermalExpansionMatproZircaloy.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef thermalExpansionMatproZircaloy_H
#define thermalExpansionMatproZircaloy_H

#include "thermalExpansionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class thermalExpansionMatproZircaloy Declaration
\*---------------------------------------------------------------------------*/

class thermalExpansionMatproZircaloy
:
    public thermalExpansionModel
{
    // Private data 

        //- Parameters for the Matpro model
        scalar par1;
        scalar par2;
        scalar par3;
        scalar par4;
        scalar par5;
        scalar par6;
        scalar par7;
        
    // Private Member Functions

        //- Disallow default bitwise copy construct
        thermalExpansionMatproZircaloy(const thermalExpansionMatproZircaloy&);

        //- Disallow default bitwise assignment
        void operator=(const thermalExpansionMatproZircaloy&);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("ZircaloyMatpro");

    // Constructors

        //- Construct from dictionary
        thermalExpansionMatproZircaloy
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dictionary& baseDict, 
            const word defaultModel   
        );

    //- Destructor
    virtual ~thermalExpansionMatproZircaloy();


    // Member Functions
    
    //- Update thermalExpansion  
    virtual void correct(symmTensorField& sf, const scalarField& T, const labelList& addr) const;

    //- Calculate alpha*T
    const symmTensor setAlphaT(const scalar& Ti) const;    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
