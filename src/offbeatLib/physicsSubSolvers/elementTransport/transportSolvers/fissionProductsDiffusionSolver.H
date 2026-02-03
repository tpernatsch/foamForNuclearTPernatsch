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
    Foam::fissionProductsDiffusionSolver

Description
    Solver for Fission products transport.

    Refs:
    Jiang W, Hales J D, Spencer B W, et al. TRISO particle fuel performance and
    failure analysis with Bison[J]. Journal of Nuclear Materials, 2021, 548: 152795.


Usage
    In solverDict:
    \verbatim
    elementTransportOptions
    {
        solvers
        (
            FpDiffusion
        );
        FpDiffusionOptions
        {
          outerPatches   outer;
        }
    }

    \endverbatim

SourceFiles
    fissionProductsDiffusionSolver.C

\mainauthor
    F. Xiang - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
    Laboratory for Reactor Physics and Systems Behaviour) & XJTU (Xi'an JiaoTong
    University, China, Nuclear Thermal-hydraulic Laboratory)

\contribution
    A. Scolaro - EPFL
    E. Brunetto - EPFL

\date
    June 2023

\*---------------------------------------------------------------------------*/

#ifndef fissionProductsDiffusionSolver_H
#define fissionProductsDiffusionSolver_H

#include "IOdictionary.H"
#include "transportSolver.H"
#include "materialModel.H"
#include "fuelMaterial.H"
#include "UO2.H"
#include "diffCoefModel.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class fissionProductsDiffusionSolver Declaration
\*---------------------------------------------------------------------------*/

class fissionProductsDiffusionSolver
:
public transportSolver
//public regIOobject
{
    // Private data

        //- Reference to temperature field
        const volScalarField& T_;

        //- Reference to fast neutron fluence field
        const volScalarField& phi_;

        //- Reference to burnup field MWd/tUO2
        const volScalarField& Bu_;

        //- Reference to volumetric heat source field
        const volScalarField& Q_;

        //- Concentration of fission products
        List<autoPtr<volScalarField>> Fp_;

        //- Thermal diffusion coefficient of different FPs
        List<autoPtr<volScalarField>> diff_;

        //- Source term for different Fps;
        List<autoPtr<volScalarField>> s_;

        List<autoPtr<diffCoefModel>> diffModels_;

        //- The flux*surface of fission products going out of fuel
        scalarList releasedInc_;

        //- The amount of fission products that have released (previous timestep)
        scalarList releasedPrev_;

        //- The amount of fission products that have released
        scalarList released_;

        scalar currentTime_;

        //- Initial amount of fission products
        scalarList totalInital_;

        //- The amount of FPs which have been produced
        scalarList totalProd_;
        scalarList totalProdPre_;

        //- The amount of FPs which have been produced
        scalarList Fpleft_;

        scalar enrichment_;


        //- Name of the boundary the Fp will escape from
        word patchName_;

        //Fission products name
        wordList FpNames_;

        //The summary name for all the fission products
        word FpName_;

    // Private Member Functions

        //- Initialize Fission products concentration
        virtual void initializeFpConcentration(const dictionary& FpOptDict);

        virtual void updateSourceTerm(const labelList& addr, const scalarField& Bu, const scalarField& f, const scalar& e);

        //- Update diffusion coefficient
        virtual void updateCoefficients();

        //- Output informations about FP
        virtual void outputInfo();

        //- Disallow default bitwise copy construct
        fissionProductsDiffusionSolver(const fissionProductsDiffusionSolver&);

        //- Disallow default bitwise assignment
        void operator=(const fissionProductsDiffusionSolver&);


public:

    //- Runtime type information
    TypeName("FpDiffusion");

    // Constructors

        //- Construct from mesh, materials and dictionary
        fissionProductsDiffusionSolver
        (
            fvMesh& mesh,
            const materials& mat,
            const dictionary& elementTransportDict,
            const word& solverName
        );

    //- Destructor
    ~fissionProductsDiffusionSolver()
    {}

    // Access Functions

        // Return summary name of fission products
        virtual const word& fieldName()
        {
            return FpName_;
        };

    // Member Functions

        //- Correct
        virtual void correctFp(volScalarField& Fp, const volScalarField diffFp, const volScalarField s);
        virtual void correct();
        //- Return true if converged
        //virtual bool converged( const word& fieldName );

        //- Next time step according to solver criterion
        virtual scalar nextDeltaT()
        {
            //Not Implemented
            return GREAT;
        };

        virtual void updateMesh(){};

        //virtual bool writeData(Ostream& os) const;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
