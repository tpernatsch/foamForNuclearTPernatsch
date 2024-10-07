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
    Foam::thermalExpansionGehr1515Ti

Description
    Class modelling the thermal expansion of 15-15 Ti cladding through the
    correlation proposed by Gehr (1973).

SourceFiles
    thermalExpansionGehr1515Ti.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef thermalExpansionGehr1515Ti_H
#define thermalExpansionGehr1515Ti_H

#include "thermalExpansionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class thermalExpansionGehr1515Ti Declaration
\*---------------------------------------------------------------------------*/

class thermalExpansionGehr1515Ti
:
    public thermalExpansionModel
{
    // Private data 

        //- Parameters for the Matpro model
        scalar par1;
        scalar par2;
        scalar par3;

        //- Perturbation parameter
        scalar perturb;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        thermalExpansionGehr1515Ti(const thermalExpansionGehr1515Ti&);

        //- Disallow default bitwise assignment
        void operator=(const thermalExpansionGehr1515Ti&);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("1515TiGehr");

    // Constructors

        //- Construct from dictionary
        thermalExpansionGehr1515Ti
        (
            const fvMesh& mesh,
            const dictionary& dict, 
            const word defaultModel   
        );

    //- Destructor
    virtual ~thermalExpansionGehr1515Ti();


    // Member Functions
    
    //- Update thermalExpansion  
    virtual void correct(symmTensorField& sf, const scalarField& T, const labelList& addr) const;   
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
