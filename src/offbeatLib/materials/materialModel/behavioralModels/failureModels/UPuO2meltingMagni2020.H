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
    Foam::UPuO2meltingMagni2020

Description
    Model to check failure of UPuO2 based on a limit melting temperature from
    "Modelling and assessment of thermal conductivity and melting
    behaviour of MOX fuel for fast reactor applications - A.Magni et al.".
    (https://doi.org/10.1016/j.jnucmat.2020.152410).

Usage


SourceFiles
    UPuO2meltingMagni2020.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\contribution
    A. Scolaro, C. Fiorina - EPFL

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef UPuO2meltingMagni2020_H
#define UPuO2meltingMagni2020_H

#include "failureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class UPuO2meltingMagni2020 Declaration
\*---------------------------------------------------------------------------*/

class UPuO2meltingMagni2020
:
    public failureModel
{
    // Private data

        //- Name used for the burnup field
        const word burnupName_;
    
        //- Reference to Burnup field
        const volScalarField* Bu_;

        //- Melting temperature for fresh UO2
        const scalar TmUO2_;

        //- OxygenMetal (O/M) ratio
        const volScalarField& OM_;

        //- Ref to number density of U (1/m3)
        const volScalarField& N_U_;

        //- Ref to number density of Pu (1/m3)
        const volScalarField& N_Pu_;

        //- Ref to number density of Am (1/m3)
        const volScalarField& N_Am_;

        //- Ref to number density of Np (1/m3)
        const volScalarField& N_Np_;

        // Correlation parameters

        // [ K/at.% ]
        const scalar gammaPu_;
        
        // [ K ]
        const scalar gammaX_;
        
        // [ K ]
        const scalar TmInf_;

        // [ GWd/tHM ]
        const scalar delta_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        UPuO2meltingMagni2020(const UPuO2meltingMagni2020&);

        //- Disallow default bitwise assignment
        void operator=(const UPuO2meltingMagni2020&);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("UPuO2Magni2020");


    // Declare run-time constructor selection table


    // Constructors

        //- Construct from dictionary
        UPuO2meltingMagni2020
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dictionary& baseDict
        );


    // Selectors


    //- Destructor
    virtual ~UPuO2meltingMagni2020();


    // Member Functions
    
        //- Check failure occurrence
        virtual void checkFailure(const labelList& addr);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
