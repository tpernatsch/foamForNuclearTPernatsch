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
    Foam::conductivityUPuO2Magni2020

Description
    Model for thermal conductivity of high-Pu content MOX. The used correlation 
    is taken from "Modelling and assessment of thermal conductivity and melting
    behaviour of MOX fuel for fast reactor applications - A.Magni et al.".
    (https://doi.org/10.1016/j.jnucmat.2020.152410)
    The ranges of applicability of the correlation are the following:
    - Temperature : [500,2700] K;
    - Deviation from stoichiometry : [0,04];
    - Plutonium content : [0,45] at.%;
    - Porosity p: [0,7] %;
    - Burnup : [0,130] GWd/tHM.

SourceFiles
    conductivityUPuO2Magni2020.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef conductivityUPuO2Magni2020_H
#define conductivityUPuO2Magni2020_H

#include "conductivityModel.H"
#include "materialModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class conductivityUPuO2Magni2020 Declaration
\*---------------------------------------------------------------------------*/

class conductivityUPuO2Magni2020
:
    public conductivityModel
{
    // Private data

        //- Name used for the burnup field
        const word burnupName_;
        
        //- Reference to Burnup field
        const volScalarField* Bu_;

        //- Asymptotic lower limit for kappa degradation with Bu [W/m K]
        scalar kappaInf_;

        //- Correlation coefficient phi [GWd/tHM]
        scalar phi_;

        //- O/M oxygen to metal ratio
        scalar OM_;

        //- Deviation from stoichiometry (2-O/M)
        scalar x_;

        //- Plutonium atomic concentration in U+Pu [at.%] 
        scalar Pu_;

        //- Weight concentration of U over U+Pu [%w]
        scalar ratioUMetal_;

        //- Weight concentration of U over U+Pu [%w]
        scalar ratioPuMetal_;

        //- U composition - weight fractions
        //- namely: (U234, U235, U236, U238)
        const scalarList isotopesU_;

        //- Pu composition - weight fractions
        //- namely: (Pu238, Pu239, Pu240, Pu241, Pu242)
        const scalarList isotopesPu_;

        // TODO (E.B.) : Do we actually need  to implement to model the porosity
        //               evolution in time? From TRANSURANUS it looks like they use 
        //               the initial porosity value for this correlation (pag 720 of
        //               TRANSURANUS manual V1M1J17) .
        //- Fuel porosity
        scalar porosity_;
        
        //- Correlation constants:
        scalar A0_;  // [m K/W]
        scalar Ax_;  // [m K/W]
        scalar APu_; // [m K/W]
        scalar B0_;  // [m/W]
        scalar BPu_; // [K/W]
        scalar D_;   // [W K/m]
        scalar E_;   // [K]

        //- Perturbation factor
        scalar perturb;
    
    // Private Member Functions

        //- Disallow default bitwise copy construct
        conductivityUPuO2Magni2020(const conductivityUPuO2Magni2020&);

        //- Disallow default bitwise assignment
        void operator=(const conductivityUPuO2Magni2020&);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("UPuO2Magni2020");

    // Constructors

        //- Construct from dictionary
        conductivityUPuO2Magni2020
        (
            const fvMesh& mesh,
            const dictionary& dict, 
            const word defaultModel   
        );
    
    //- Destructor
    virtual ~conductivityUPuO2Magni2020();


    // Member Functions
    
        //- Update conductivity  
        virtual void correct(scalarField& sf, const scalarField& T, const labelList& addr);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
