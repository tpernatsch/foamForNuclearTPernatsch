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
    Foam::MonolithicSiCCreepModel

Description
    Compute creep strain rate in the following way:
    creep_rate = K*sigmaEff*fast_neutron_flux;

    K can be provided by users, otherwised it will be calculated by the correlation.

Usage
                \verbatim
                rheology byMaterial;

                ...

                // List of materials, one per cellZone.
                materials
                {
                    ...

                    "SiC"
                    {
                      material SiC;
                      ...
                      rheologyModel  misesPlasticCreep;
                      rheologyModelOptions
                      {
                          plasticStrainVsYieldStress table
                          (
                              (0    1e60)
                          );
                          relax 1;
                          creepModel  SiCMonolithic;
                          fluxConversionFactor       1.0;
                      }

                    ...
                }
                \endverbatim

SourceFiles
    MonolithicSiCCreepModel.C

\mainauthor
    F. Xiang - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
    Laboratory for Reactor Physics and Systems Behaviour) & XJTU (Xi'an JiaoTong
    University, China, Nuclear Thermal-hydraulic Laboratory)

\contribution
    A. Scolaro - EPFL
    E. Brunetto - EPFL

\date
    December 2022

\*---------------------------------------------------------------------------*/

#ifndef ParfumeBufferCreepModel_H
#define ParfumeBufferCreepModel_H

#include "creepModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class thermoMechanicsSolver Declaration
\*---------------------------------------------------------------------------*/

class MonolithicSiCCreepModel
:
    public creepModel
{
    // Private Member Functions

        scalar timeIndex_;

        scalar relax_;

        //- Name used for the fastFlux field
        const word fastFluxName_;

        //- Creep coefficient.
        const scalar K_;

        //-Poisson’s ratio for creep
        const scalar vc_;

        //- Flux conversion factor. Normally the fast neutron is defined as the neutron
        //- with energy(E) more than 0.1 MeV. However in this model, the neutron (E > 0.18 MeV)
        //- flux is used. So the conversion between flux(>0.1 MeV) to flux(>0.18 MeV) is needed.
        const scalar fcf_;

        //- Model variables
        const scalar par1;
        const scalar par2;
        const scalar par3;
        const scalar par4;
        const scalar par5;

        //- Disallow default bitwise copy construct
        MonolithicSiCCreepModel(const MonolithicSiCCreepModel&);

        //- Disallow default bitwise assignment
        void operator=(const MonolithicSiCCreepModel&);

public:

    //- Runtime type information
    TypeName("SiCMonolithic");

    // Constructors

        //- Construct from mesh and thermo-mechanical properties
        MonolithicSiCCreepModel
        (
            const fvMesh& mesh,
            const dictionary& lawDict
        );

    //- Destructor
    ~MonolithicSiCCreepModel();


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
