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
    Foam::MatproCreepModel

Description
    Matpro creep model, derived from BISON manual and SCDAP/RELAP5/MOD 3.3
    Code Manual.

    \htmlinclude creepModel_MatproCreepModel.html

SourceFiles
    MatproCreepModel.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)\n

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef MatproCreepModel_H
#define MatproCreepModel_H

#include "creepModel.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class "MatproCreepModel" Declaration
\*---------------------------------------------------------------------------*/

class MatproCreepModel
:
    public creepModel
{
    // Private Member Functions    
    
        //- Names for density, grain radius, heatSource and oxygen metal ratio fields
        const word heatSourceName_;
        const word densityName_;        
        const word grainRadiusName_;        
        const word oxygenMetalRatioName_;

        //- Apply correction from Sakai (present also in BISON)
        Switch SakaiCorrection_;

        //- Model variables
        const scalar A1;
        const scalar A2;
        const scalar A3;
        const scalar A4;
        const scalar A6;
        const scalar A7;

        //- Universal Gas constant
        const scalar R;

        scalar timeIndex_;

        scalar relax_;
        
        //- Disallow default bitwise copy construct
        MatproCreepModel(const MatproCreepModel&);

        //- Disallow default bitwise assignment
        void operator=(const MatproCreepModel&);

public:
 
    //- Runtime type information
    TypeName("MatproCreepModel");


    // Constructors

        //- Construct from mesh and thermo-mechanical properties
        MatproCreepModel
        (
            const fvMesh& mesh,
            const dictionary& lawDict
        );


    //- Destructor
    ~MatproCreepModel();


    // Member Functions
    
        //- Update mechanical properties according to rheology model        
        virtual void correctCreep
        (
            volSymmTensorField& epsilonEl,
            const labelList& addr
        );         
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
