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
    Foam::conductivityUPuO2Philipponneau

Description
    Model for thermal conductivity of MOX taken from "Thermal conductivity of 
    (U,Pu)O2−x mixed oxide fuel, Y. Philiponneau, 1992"
    (https://doi.org/10.1016/0022-3115(92)90470-6).

    Note: in this correlation the MOX conductivity is independent of the Pu 
    content.

    Note 2: this correlation is valid for temperatures 500 K < T < melting T.

    Note 3: as suggested in (https://doi.org/10.1016/j.pnucene.2017.03.016) the 
    conversion factor from Bu[GWd/tHM] to FIMA[%] in MOX is taken equal to 9.5 .

SourceFiles
    conductivityUPuO2Philipponneau.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    August 2022

\*---------------------------------------------------------------------------*/

#ifndef conductivityUPuO2Philipponneau_H
#define conductivityUPuO2Philipponneau_H

#include "conductivityModel.H"
#include "materialModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class conductivityUPuO2Philipponneau Declaration
\*---------------------------------------------------------------------------*/

class conductivityUPuO2Philipponneau
:
    public conductivityModel
{
    // Private data

        //- Name used for the burnup field
        const word burnupName_;
        
        //- Reference to Burnup field
        const volScalarField* Bu_;

        //- B [ m/W ]
        scalar B_;

        //- C [ W/(m K^4) ]
        scalar C_;

        //- OxygenMetal (O/M) ratio
        const volScalarField& OM_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        conductivityUPuO2Philipponneau(const conductivityUPuO2Philipponneau&);

        //- Disallow default bitwise assignment
        void operator=(const conductivityUPuO2Philipponneau&);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("UPuO2Philipponeau");

    // Constructors

        //- Construct from dictionary
        conductivityUPuO2Philipponneau
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dictionary& baseDict,
            const word defaultModel         
        );
    
    //- Destructor
    virtual ~conductivityUPuO2Philipponneau();


    // Member Functions
    
        //- Update conductivity  
        virtual void correct(scalarField& sf, const scalarField& T, const labelList& addr);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
