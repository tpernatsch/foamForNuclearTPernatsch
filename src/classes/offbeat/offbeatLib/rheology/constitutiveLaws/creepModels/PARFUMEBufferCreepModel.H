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
    Foam::PARFUMEBufferCreepModel

Description
    Compute creep strain rate in the following way:
    creep_rate = K*[s1 - vc * (s2 + s3)]*fast_neutron_flux;
    s1, s2, s3 - Principal stress;
    vc - Poisson ratio for creep.

    Usage
                \verbatim
                rheology byMaterial;

                ...

                // List of materials, one per cellZone.
                materials
                {
                    ...

                    "IPyC|OPyC"
                    {
                      material PyC;
                      ...
                      rheologyModel misesPlasticCreep;

                      rheologyModelOptions
                      {
                          plasticStrainVsYieldStress table
                          (
                              (0    1e60)
                          );
                           relax 1.0;
                          creepModel PARFUMUEPyCCreepModel;
                          fluxConversionFactor       0.91;
                      }

                    ...
                }
                \endverbatim

SourceFiles
    PARFUMEBufferCreepModel.C

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

#ifndef PARFUMEBufferCreepModel_H
#define PARFUMEBufferCreepModel_H

#include "creepModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class thermoMechanicsSolver Declaration
\*---------------------------------------------------------------------------*/

class PARFUMEBufferCreepModel
:
    public creepModel
{
    // Private Member Functions

        scalar timeIndex_;

        scalar relax_;

        //- Name for density field
        const word densityName_;
        //- Name used for the fastFlux field
        const word fastFluxName_;

        //- Creep amplification coefficient, the default value is 2.0 to better fit
        //- experimental data.
        const scalar C_amp;

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
        PARFUMEBufferCreepModel(const PARFUMEBufferCreepModel&);

        //- Disallow default bitwise assignment
        void operator=(const PARFUMEBufferCreepModel&);

public:

    //- Runtime type information
    TypeName("PARFUMUEBufferCreepModel");

    // Constructors

        //- Construct from mesh and thermo-mechanical properties
        PARFUMEBufferCreepModel
        (
            const fvMesh& mesh,
            const dictionary& lawDict
        );

    //- Destructor
    ~PARFUMEBufferCreepModel();


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
