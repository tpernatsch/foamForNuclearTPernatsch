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
    Foam::constantCreepPrincipalStress

Description
    Compute creep strain rate based on the eigen values of the stress tensor.
    This class is used to verify the OFFBEAT results with IAEA CPR-6 benchmark for TRISO fuel applications.
    
    \htmlinclude creepModel_constantCreepPrincipalStress.html

SourceFiles
    constantCreepPrincipalStress.C

\mainauthor
    F. Xiang - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
    Laboratory for Reactor Physics and Systems Behaviour) & XJTU (Xi'an JiaoTong
    University, China, Nuclear Thermal-hydraulic Laboratory)

\contribution
    A. Scolaro - EPFL
    E. Brunetto - EPFL

\date
    January 2023

\*---------------------------------------------------------------------------*/

#ifndef PARFUMEBufferCreepModel_H
#define PARFUMEBufferCreepModel_H

#include "creepModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class thermoMechanicsSolver Declaration
\*---------------------------------------------------------------------------*/

class constantCreepPrincipalStress
:
    public creepModel
{
    // Private Member Functions

        scalar timeIndex_;

        scalar relax_;

        //- Name used for the fastFlux field
        const word fastFluxName_;

        //-Poisson’s ratio for creep
        const scalar vc_;

        //- Flux conversion factor. Normally the fast neutron is defined as the neutron
        //- with energy(E) more than 0.1 MeV. However in this model, the neutron (E > 0.18 MeV)
        //- flux is used. So the conversion between flux(>0.1 MeV) to flux(>0.18 MeV) is needed.
        const scalar fcf_;

        //- Creep coefficient. Provided by users in unit of /(MPa 10^25n/m²)
        const scalar K_;


        //- Disallow default bitwise copy construct
        constantCreepPrincipalStress(const constantCreepPrincipalStress&);

        //- Disallow default bitwise assignment
        void operator=(const constantCreepPrincipalStress&);

public:

    //- Runtime type information
    TypeName("constantCreepPrincipalStress");

    // Constructors

        //- Construct from mesh and thermo-mechanical properties
        constantCreepPrincipalStress
        (
            const fvMesh& mesh,
            const dictionary& lawDict
        );

    //- Destructor
    ~constantCreepPrincipalStress();


    // Member Functions


        //- Update mechanical properties according to rheology model
        virtual void correctCreep
        (
            volSymmTensorField& epsilonEl,
            const labelList& addr
        );

        //- Compute the principal stress. If they can not be solve, the value 0.0
        //- will be used and solved = false. stressPrinciple1 > stressPrinciple2 > stressPrinciple3
        void computePrincipalStress
        (
            symmTensor s,
            scalar& stressPrinciple1,
            scalar& stressPrinciple2,
            scalar& stressPrinciple3,
            bool& solved
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
