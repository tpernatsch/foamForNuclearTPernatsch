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
    Foam::conductivityUPuO2MagniMAMOX

Description
    Model for thermal conductivity of high-Pu content MOX. The used correlation 
    is taken from "Modelling of thermal conductivity and melting behaviour of 
    minor actinide-MOX fuels and assessment against experimental and molecular 
    dynamics data - A.Magni et al.".
    (https://doi.org/10.1016/j.jnucmat.2021.153312)
    
SourceFiles
    conductivityUPuO2MagniMAMOX.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef conductivityUPuO2MagniMAMOX_H
#define conductivityUPuO2MagniMAMOX_H

#include "conductivityModel.H"
#include "materialModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class conductivityUPuO2MagniMAMOX Declaration
\*---------------------------------------------------------------------------*/

class conductivityUPuO2MagniMAMOX
:
    public conductivityModel
{
    // Private data

        //- Name used for the burnup field
        const word burnupName_;
        
        //- Reference to Burnup field
        const volScalarField* Bu_;

        //- OxygenMetal (O/M) ratio
        const volScalarField& OM_;

        //- Ref to number density of U (1/m3)
        const volScalarField& N_U_;

        //- Ref to number density of Pu (1/m3)
        const volScalarField& N_Pu_;

        //- Ref to number density of Am (1/m3)
        volScalarField* N_Am_;

        //- Ref to number density of Np (1/m3)
        volScalarField* N_Np_;

        //- Fuel porosity
        const volScalarField& porosity_;

        //- Asymptotic lower limit for kappa degradation with Bu [W/m K]
        scalar kappaInf_;

        //- Correlation coefficient phi [GWd/tHM]
        scalar phi_;

        //- Correlation constants:
        scalar A0_;  // [m K/W]
        scalar Ax_;  // [m K/W]
        scalar APu_; // [m K/W]
        scalar AAm_; // [m K/W]
        scalar ANp_; // [m K/W]
        scalar B0_;  // [m/W]
        scalar BPu_; // [K/W]
        scalar BAm_; // [K/W]
        scalar BNp_; // [K/W]
        scalar D_;   // [W K/m]
        scalar E_;   // [K]
    
    // Private Member Functions

        //- Disallow default bitwise copy construct
        conductivityUPuO2MagniMAMOX(const conductivityUPuO2MagniMAMOX&);

        //- Disallow default bitwise assignment
        void operator=(const conductivityUPuO2MagniMAMOX&);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("MaUPuO2Magni");

    // Constructors

        //- Construct from dictionary
        conductivityUPuO2MagniMAMOX
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dictionary& baseDict,
            const word defaultModel         
        );
    
    //- Destructor
    virtual ~conductivityUPuO2MagniMAMOX();


    // Member Functions
    
        //- Update conductivity  
        virtual void correct(scalarField& sf, const scalarField& T, const labelList& addr);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
