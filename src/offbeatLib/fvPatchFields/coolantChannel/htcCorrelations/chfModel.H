/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
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
    along with OpenFOAM. If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::chfModel

Description
    Critical Heat Flux (CHF) model.

    Supported models:
    - EPRI           : EPRI (1985) look-up-based correlation for PWR rods
    - GroeneveldTable: Groenveld et al. look-up table (read from chfTableDict)
    - Zuber          : Zuber (1958) pool boiling limit
    - ModifiedZuber  : Zuber with subcooling correction (Frederking)
    - GE             : GE (1973) simple quality-based limit

    The GroeneveldTable model requires a constant/chfTableDict file.
    Initialise the table data by calling initializeTable() after construction.

    T_CHF is NOT computed here; it is the responsibility of the BC
    (via solvePreCHFTemperature or the fast approximation).

SourceFiles
    chfModel.C

\*---------------------------------------------------------------------------*/

#ifndef chfModel_H
#define chfModel_H

#include "scalar.H"
#include "word.H"
#include "label.H"
#include "scalarField.H"
#include "FieldField.H"
#include "PtrList.H"
#include "autoPtr.H"
#include "InterpolateTables.H"
#include "PtrListScalarFieldFieldINew.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class dictionary;
class fvMesh;

/*---------------------------------------------------------------------------*\
                         Class chfModel Declaration
\*---------------------------------------------------------------------------*/

class chfModel
{
public:

    enum class Model
    {
        EPRI,
        GroeneveldTable,
        Zuber,
        ModifiedZuber,
        GE
    };


    // Constructors

        chfModel(const word& name);

        //- Copy constructor: rebuilds the interpolation table from copied data
        chfModel(const chfModel& other);


    // Member Functions

        //- Pre-compute non-local patch quantities needed by some models.
        //  Must be called once per timestep before the per-face compute loop.
        //  Non-op for all models except EPRI.
        void precompute
        (
            const scalarField& zetas,           // face axial positions   [m]
            const scalarField& faceHeights,     // face height weights    [m]
            const scalarField& heatFlux,        // wall heat flux         [W/m²]
            const scalarField& massFlowRate,    // per-face mass flow rate [kg/s]
            const scalarField& coolantPressure, // per-face pressure      [Pa]
            scalar inletTemperature,            // inlet coolant temp     [K]
            scalar flowArea,                    // flow area              [m²]
            scalar flowDirection                // +1 or -1
        );

        //- Return q''_CHF [W/m²] for face faceI.
        //  Call precompute() before the face loop.
        scalar compute
        (
            scalar p,       // pressure                              [Pa]
            scalar T0,      // bulk coolant temperature (ModifiedZuber) [K]
            scalar Q,       // thermodynamic quality                 [-]
            scalar G,       // mass flux                             [kg/(m²·s)]
            scalar DH,      // hydraulic diameter                    [m]
            scalar qWall,   // current wall heat flux                [W/m²]
            label  faceI    // face index (EPRI uses FAP_[faceI])
        ) const;

        //- Initialise the Groenveld look-up table (GroeneveldTable model only).
        //  Reads chfTableDict from constant/ and interpolation options from bcDict.
        void initializeTable(const dictionary& bcDict, const fvMesh& mesh);

        Model model() const { return model_; }


private:

    Model model_;

    // EPRI pre-computed non-local quantities (set by precompute())

        //- Inlet thermodynamic quality (EPRI)
        scalar Qin_;

        //- Per-face axial flux shape factor (EPRI)
        scalarField FAP_;

    // Private helper methods (EPRI)

        //- Compute inlet quality from inlet temperature and pressure
        static scalar computeInletQuality
        (
            scalar inletTemperature,
            scalar pInlet
        );

        //- Compute per-face axial flux shape factors
        static scalarField computeShapeFactors
        (
            const scalarField& zetas,
            const scalarField& faceHeights,
            const scalarField& heatFlux,
            const scalarField& massFlowRate,
            scalar flowArea,
            scalar flowDirection
        );

    // Groenveld table data (only used when model_ == GroeneveldTable)

        PtrList<FieldField<Field, scalar>> chfData_;

        scalar chfPitchToDiameter_;

        scalarField chfPressureValues_;
        scalarField chfMassFluxValues_;
        scalarField chfQualityValues_;

        interpolateTableBase::interpolationMethod chfPressureMethod_;
        interpolateTableBase::interpolationMethod chfMassFluxMethod_;
        interpolateTableBase::interpolationMethod chfQualityMethod_;

        autoPtr<scalarFieldFieldInterpolateTable> chfPressureTable_;

    // Unit conversion constants

        static const scalar PsiToPa_;
        static const scalar BtuToJ_;
        static const scalar FeetToMeters_;
        static const scalar LbToKg_;
        static const scalar HrToSec_;

        static const scalar g_;

    static Model parse(const word& name);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
