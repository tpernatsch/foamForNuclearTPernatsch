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
    Foam::thermalExpansionLemehovUPuO2

Description
    Class modelling thermal expansion of (U,Pu)O2 MOX fuel from Lemehov(2020).
    This correlation is valid for:
    - O/M ratios between 1.94 and 2;
    - Pu/HM ratios lower than 60 %.

    Source : https://re.public.polimi.it/retrieve/handle/11311/1172415/619939/Technical%20Report_INSPYRE_WP7-D7.2_%282020%29.pdf

SourceFiles
    thermalExpansionLemehovUPuO2.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef thermalExpansionLemehovUPuO2_H
#define thermalExpansionLemehovUPuO2_H

#include "thermalExpansionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class thermalExpansionLemehovUPuO2 Declaration
\*---------------------------------------------------------------------------*/

class thermalExpansionLemehovUPuO2
:
    public thermalExpansionModel
{
    // Private data

        //- Parameters for the Matpro model
        scalar b0_;
        scalar b1_;
        scalar b2_;
        scalar b3_;
        scalar by_;

        //- O/M ratio
        const volScalarField& OM_;
        
        //- Ref to number density of U (1/m3)
        const volScalarField& N_U_;

        //- Ref to number density of Pu (1/m3)
        const volScalarField& N_Pu_;

        //- Ref to number density of Am (1/m3)
        volScalarField* N_Am_;

        //- Ref to number density of Np (1/m3)
        volScalarField* N_Np_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        thermalExpansionLemehovUPuO2(const thermalExpansionLemehovUPuO2&);

        //- Disallow default bitwise assignment
        void operator=(const thermalExpansionLemehovUPuO2&);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("UPuO2Lemehov");

    // Constructors

        //- Construct from dictionary
        thermalExpansionLemehovUPuO2
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dictionary& baseDict, 
            const word defaultModel   
        );

    //- Destructor
    virtual ~thermalExpansionLemehovUPuO2();


    // Member Functions
    
    //- Update thermalExpansion  
    virtual void correct(symmTensorField& sf, const scalarField& T, const labelList& addr) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
