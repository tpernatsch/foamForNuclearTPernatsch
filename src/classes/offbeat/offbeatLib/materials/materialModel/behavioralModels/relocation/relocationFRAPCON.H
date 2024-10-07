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
    Foam::relocationFRAPCON

Description
    Class modelling the relocation phoenomenon derived from FRAPCON.

Usage
    In solverDict file:
    \verbatim
    materials
    {
        fuel
        {
            material UO2;

            ...

            relocationModel   UO2FRAPCON;
            GapCold           1e-5;
            DiamCold          1e-2;
            outerPatch        fuelOuterSurface;
        }

        ...
    }
    \endverbatim

SourceFiles
    relocationFRAPCON.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef relocationFRAPCON_H
#define relocationFRAPCON_H

#include "relocationModel.H"
#include "physicoChemicalConstants.H"
#include "sliceMapper.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class relocationFRAPCON Declaration
\*---------------------------------------------------------------------------*/

class relocationFRAPCON
:
    public relocationModel
{
    // Private data
        
        //- Name used for the burnup field
        const word burnupName_;
        
        //- Name used for the heat Source field
        const word heatSourceName_;
        
        //- Name used for the gapWidth field
        const word gapWidthName_;

        //- Fuel outersurface patchID to estimate average gapwidth
        label outerPatchID_;
        
        //- Reference to gapWidth field
        const volScalarField* gapWidth_;
        
        //- Reference to slice mapper
        const sliceMapper* mapper_;
        
        //- Cold Gap Reference Thickness [m]
        scalar GapCold_;

        //- Cold Pellet Diamenter [m]
        scalar DiamCold_;
        
        //- Maximum fraction of relocation recovered
        scalar recFrac_;

        //- Perturbation parameter
        scalar perturb;

        //- Decide whether to use modified version of FRAPCON model 
        //- used since FRAPCON 3.5 (on by default)
        Switch modifiedRelocationModel_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        relocationFRAPCON(const relocationFRAPCON&);

        //- Disallow default bitwise assignment
        void operator=(const relocationFRAPCON&);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("UO2FRAPCON");


    // Declare run-time constructor selection table


    // Constructors

        //- Construct from dictionary
        relocationFRAPCON
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

    // Selectors
        

    //- Destructor
    virtual ~relocationFRAPCON();


    // Member Functions
    
        //- Update relocation strain  
        virtual void correct(const scalarField& T, const labelList& addr);
        
        //- Set relocation strain
        Foam::scalar setRelocation(double Q, double Gap, double D, double Bu);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
