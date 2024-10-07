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
    Foam::thermalExpansionMatproUPuO2

Description
    Class modelling thermal expansion of (U,Pu)O2 MOX fuel from MATPROv11.

SourceFiles
    thermalExpansionMatproUPuO2.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef thermalExpansionMatproUPuO2_H
#define thermalExpansionMatproUPuO2_H

#include "thermalExpansionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class thermalExpansionMatproUPuO2 Declaration
\*---------------------------------------------------------------------------*/

class thermalExpansionMatproUPuO2
:
    public thermalExpansionModel
{
    // Private data

        //- Parameters for the MATPRO model
        scalar par1_;
        scalar par2_;
        scalar par3_;
        scalar par4_;

        //- Perturbation parameter
        scalar perturb;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        thermalExpansionMatproUPuO2(const thermalExpansionMatproUPuO2&);

        //- Disallow default bitwise assignment
        void operator=(const thermalExpansionMatproUPuO2&);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("UPuO2MATPRO");

    // Constructors

        //- Construct from dictionary
        thermalExpansionMatproUPuO2
        (
            const fvMesh& mesh,
            const dictionary& dict, 
            const word defaultModel   
        );

    //- Destructor
    virtual ~thermalExpansionMatproUPuO2();


    // Member Functions
    
    //- Update thermalExpansion  
    virtual void correct(symmTensorField& sf, const scalarField& T, const labelList& addr) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
