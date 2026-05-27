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
    Solver for fission product (FP) transport by diffusion.

    Solves a transient diffusion equation for each configured FP species:
    \f[ \frac{\partial C}{\partial t} - \nabla \cdot (D \nabla C) = S \f]

    Diffusion coefficients are provided by material-specific \c diffCoefModel
    subclasses.  Source terms (fission yields) are provided by per-species
    \c fpYieldModel subclasses.

    Refs:
    Jiang W, Hales J D, Spencer B W, et al. TRISO particle fuel performance and
    failure analysis with Bison[J]. Journal of Nuclear Materials, 2021, 548: 152795.


Usage
    In solverDict (material-first format):
    \verbatim
    elementTransportOptions
    {
        solvers ( FpDiffusion );
        FpDiffusionOptions
        {
            fissionProducts (Cs Ag);     // word list of species

            diffusion                    // zone-keyed; wordRe patterns accepted
            {
                Kernel      { type UO2Arrhenius; }
                Buffer      { type singleArrhenius; d1 1e-8; q1 0; }
                "IPyC|OPyC" { type PyCArrhenius; }
                SiC         { type SiCArrhenius; }
            }

            yield                        // zone-keyed; only fissile zones needed
            {
                Kernel { type UO2; }     // smart default: built-in Cs/Ag/Kr/Sr yields
            }

            // Per-species overrides — add a species subdict inside the zone dict:
            // diffusion
            // {
            //     Kernel
            //     {
            //         Cs { type doubleArrhenius; d1 5.6e-8; q1 209000; d2 5.2e-4; q2 362000; }
            //         Ag { type UO2Arrhenius; }
            //     }
            // }
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

#include "transportSolver.H"
#include "materialModel.H"
#include "diffCoefModel.H"
#include "fpYieldModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class fissionProductsDiffusionSolver Declaration
\*---------------------------------------------------------------------------*/

class fissionProductsDiffusionSolver
:
    public transportSolver
{
    // Private data

        //- Reference to temperature field
        const volScalarField& T_;

        //- Fast neutron fluence field (nullptr when not in registry)
        const volScalarField* phi_;

        //- Burnup field [MWd/tUO2] (nullptr when not in registry)
        const volScalarField* Bu_;

        //- Reference to volumetric heat source field
        const volScalarField& Q_;

        //- Concentration fields, one per FP species [mol/m³]
        List<autoPtr<volScalarField>> Fp_;

        //- Diffusion coefficient fields, one per FP species [m²/s]
        List<autoPtr<volScalarField>> diff_;

        //- Source term fields, one per FP species [mol/m³/s]
        List<autoPtr<volScalarField>> s_;

        //- Diffusion coefficient models [species][material zone]
        List<List<autoPtr<diffCoefModel>>> diffModels_;

        //- Yield (source) models, one per FP species
        List<autoPtr<fpYieldModel>> yieldModels_;

        //- Release rate (from mass balance) [mol/s]
        scalarList releasedInc_;

        //- Cumulative release at the previous timestep [mol]
        scalarList releasedPrev_;

        //- Cumulative release [mol]
        scalarList released_;

        //- Time index at which releasedPrev_ / totalProdPre_ were last stored
        label timeIndex_;

        //- Initial inventory [mol]
        scalarList totalInitial_;

        //- Cumulative production [mol]
        scalarList totalProd_;
        scalarList totalProdPrev_;

        //- Remaining inventory in domain [mol]
        scalarList Fpleft_;

        //- Enrichment of the UO2 kernel (atom fraction)
        scalar enrichment_;

        //- FP species names (keys of fissionProducts subdict)
        wordList FpNames_;

    // Private Member Functions

        //- Update diffusion coefficients and source terms for all species
        void updateCoefficients();

        //- Output release / production statistics to Info
        void outputInfo();

        fissionProductsDiffusionSolver(const fissionProductsDiffusionSolver&) = delete;
        void operator=(const fissionProductsDiffusionSolver&) = delete;


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
    ~fissionProductsDiffusionSolver() = default;

    // Access

        virtual const word& fieldName()
        {
            return FpNames_[0];
        }

        using transportSolver::converged;

    // Member Functions

        //- Solve diffusion for a single FP species
        void correctFp
        (
            volScalarField& Fp,
            const volScalarField& diffFp,
            const volScalarField& s
        );

        //- Correct all FP species
        virtual void correct();

        //- Return true if all FP species are converged
        virtual bool converged()
        {
            forAll(FpNames_, i)
            {
                if (!transportSolver::converged("N_" + FpNames_[i]))
                    return false;
            }
            return true;
        }

        virtual scalar nextDeltaT()
        {
            return GREAT;
        }

        virtual void updateMesh()
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
