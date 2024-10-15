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
    Foam::densificationFRAPCON

Description
    Class modelling the densification phenomenon derived from FRAPCON.

Usage
    In solverDict file:
    \verbatim
    materials
    {
        fuel
        {
            material UO2;

            ...

            densificationModel         UO2FRAPCON;
            resinteringDensityChange   0.5;

            //- Resintering temperature, 1800K by default
            //  Used only if resintering density change is 0.0
            Tsintering      1800;
        }

        ...
    }
    \endverbatim

SourceFiles
    densificationFRAPCON.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\contribution
    E. Brunetto, C. Fiorina - EPFL

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef densificationFRAPCON_H
#define densificationFRAPCON_H

#include "densificationModel.H"
#include "physicoChemicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class densificationFRAPCON Declaration
\*---------------------------------------------------------------------------*/

class densificationFRAPCON
:
    public densificationModel
{
    // Private data

        //- Name used for the burnup field
        const word burnupName_;
        
        //- Reference to Burnup field
        const volScalarField* Bu_;
        
        //- Fuel density as fraction of TD
        scalar densityFrac_;        
        
        //- Resintereing density change (fraction of %TD)
        scalar resintDensChange_;

        //- Sintering temperature
        scalar Tsinter_;

        //- Parameters for the FRAPCON model
        scalar par1;
        scalar par2;
        scalar par3;
        scalar par4;
        scalar par5;
        scalar par6;

        //- Perturbation parameter
        scalar perturb;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        densificationFRAPCON(const densificationFRAPCON&);

        //- Disallow default bitwise assignment
        void operator=(const densificationFRAPCON&);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("UO2FRAPCON");

    
    // Declare run-time constructor selection table

    
    // Constructors

        //- Construct from dictionary
        densificationFRAPCON
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

    // Selectors

        //- Destructor
        virtual ~densificationFRAPCON();


    // Member Functions
    
        //- Update heatCapacity  
        virtual void correct(const scalarField& T, const labelList& addr);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
