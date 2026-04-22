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
    Foam::coolantChannelRIAfvPatchScalarField

Description
    Coolant channel boundary condition for Reactivity Initiated Accident (RIA)
    transients.

    The model applies the Bessiron (2007) criterion for boiling crisis: the
    crisis is triggered when the vaporized water layer thickness exceeds a
    threshold, rather than when the wall heat flux exceeds a critical value.
    The bulk coolant temperature is prescribed as a time/space-dependent axial
    profile (enthalpyModel is always disabled for RIA scenarios).

    The boiling curve is divided into four regimes driven by the vaporized
    layer thickness:

    - **Pre-crisis heat conduction** (vaporized thickness below threshold,
      Tw < TCHF): transient heat conduction in stagnant water.

    - **Vaporizing** (vaporized thickness below threshold, Tw >= TCHF):
      the effective wall temperature is held at TCHF and the vaporized water
      layer thickness is integrated. When it exceeds thicknessBessironRIA, the
      thickness is capped and qCHF is frozen.

    - **Transition boiling** (vaporized thickness at threshold, Tsat < Tw):
      Bessiron exponential formula:
      qTB = qFilm + (qCHF - qFilm)*exp(-c*(Tw - TCHF)^2)

    - **Natural convection after liquid contact** (vaporized thickness at
      threshold, Tw <= Tsat): natural convection is used. The vaporized layer
      and conduction clock are reset only after the wall cools close to the
      coolant temperature.

Usage
    Example:
    \verbatim
        cladOuter
        {
            type                        coolantChannelRIA;
            value                       uniform 560;

            kappa                       k;
            hydraulicDiameter           0.01;
            flowArea                    1e-4;
            relax                       0.15;

            // Bessiron RIA parameters
            pulseStartTime               100.0;  // user time at which pulse begins
            filmBoilingRIACorrection         5.0;
            filmBoilingRIACorrectionDuration 15.0;
            leidenfrostdTBessiron            450.0;
            onbdTBessiron                    20.0;
            thicknessBessironRIA             3e-5;
            resetDeltaTBessiron              1.0;

            correlations
            {
                singlePhase     ChurchillChu;
                filmBoiling     Sakurai;
            }

            axialProfileDict
            {
                timePoints      ( 0.0 1.0 );
                axialLocations  ( 0.0 0.5 1.0 );
                T0Data          ( ( 560 560 560 )
                                  ( 565 567 570 ) );
            }

            coolantPressure     uniform 1.55e7;
            massFlowRate        uniform 0.0;
        }
    \endverbatim

SourceFiles
    coolantChannelRIAfvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef coolantChannelRIAfvPatchScalarField_H
#define coolantChannelRIAfvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

#include "sliceMapper.H"

#include "coolantChannelModel/coolantChannelModel.H"
#include "singlePhaseHTCModel.H"
#include "filmBoilingHTCModel.H"

namespace Foam
{

class coolantChannelRIAfvPatchScalarField
:
    public fvPatchField<scalar>
{
private:

    // Private data

        //- Name of thermal conductivity field in the solid region
        word kappaName_;

        //- Correlation selection dictionary
        dictionary correlations_;

        // Heat transfer regime models
        singlePhaseHTCModel singlePhase_;
        filmBoilingHTCModel filmBoiling_;

        //- Coolant channel properties and state (prescribed T0 mode only)
        coolantChannelModel coolant_;

        //- Heat transfer coefficient at patch faces
        scalarField h_;

        //- Mixed BC coefficient
        scalarField alpha_;

        //- Relaxation factor for HTC update
        scalar relax_;

        //- Auxiliary time marker used to control time-advance updates
        scalar timeCheck_;

        //- Wall heat flux
        scalarField heatFlux_;

        //- Wall temperature
        scalarField Twall_;

        //- Face heights used for axial coolant temperature interpolation
        scalarField faceHeights_;


        // --- Bessiron RIA parameters ---

        //- Film boiling HTC multiplier at the moment of boiling crisis (default 5.0).
        //  Linearly decreases to 1 over filmBoilingRIACorrectionDuration_ seconds
        //  of cumulative film boiling time.
        scalar filmBoilingRIACorrection_;

        //- Time over which filmBoilingRIACorrection_ decays from its initial
        //  value to 1 [s] (default 15.0)
        scalar filmBoilingRIACorrectionDuration_;

        //- Temperature increment above Tsat defining TLeidenfrost [K] (default 450)
        scalar leidenfrostdTBessiron_;

        //- Temperature increment above Tsat defining TCHF [K] (default 20)
        scalar onbdTBessiron_;

        //- Vaporized layer thickness that triggers boiling crisis [m] (default 3e-5)
        scalar thicknessBessironRIA_;

        //- Temperature margin above coolant temperature for post-crisis reset [K]
        scalar resetDeltaTBessiron_;

        //- User time at which the RIA pulse starts.
        //  The transient conduction clock (conductionTime_) only begins
        //  advancing from this moment; before it, single-phase NC is used.
        //  Default 0.0 (pulse starts at simulation start).
        scalar pulseStartTime_;


        // --- RIA transient state (persistent across timesteps) ---

        //- Vaporized water layer thickness [m]
        scalarField vaporizedThickness_;

        //- Vaporized water layer thickness at previous timestep [m]
        scalarField vaporizedThicknessOld_;

        //- Critical heat flux frozen at the moment of boiling crisis [W/m²]
        scalarField qCHF_;

        //- Heat flux at the Leidenfrost point [W/m²]
        scalarField qLeidenfrost_;

        //- Bessiron decay constant k [K^-2], frozen at crisis.
        //  Derived from the bare Sakurai curve (fbCorr=1) so that the
        //  minimum of the TB formula sits at TLeidenfrost in steady state.
        scalarField kBessiron_;

        //- Time spent in pre-crisis stagnant-water conduction [s]
        scalarField conductionTime_;

        //- Time spent in pre-crisis stagnant-water conduction at previous timestep [s]
        scalarField conductionTimeOld_;

        //- Cumulative time spent in post-crisis film/transition boiling [s]
        scalarField filmBoilingTime_;

        //- filmBoilingTime_ at previous timestep [s]
        scalarField filmBoilingTimeOld_;


    // Private Member Functions

        //- Compute patch face heights projected along the pin direction
        scalarField calculateFaceHeights() const;

        //- Validate geometry and required user inputs
        void validateGeometry() const;

        //- Update face heights when the mesh has moved
        void updateFaceHeightsIfNeeded();

        //- Compute solid-side conductance alphaP [W/(m²·K)]
        scalarField computeSolidSideConductance() const;

        //- Compute wall heat flux from current wall temperature estimate
        void computeHeatFlux(const scalarField& alphaP);

        //- Advance per-face vaporized thickness and flags for a new timestep
        void advanceTimeStep();

        //- Compute updated vaporized thickness for face faceI
        scalar computeVaporizedThickness
        (
            const label faceI,
            const scalar heatFlux,
            const scalar hSP,
            const scalar Tw,
            const scalar T0,
            const scalar p,
            const scalar deltaT
        ) const;

        //- Derive Bessiron decay constant k [K^-2] from the condition that
        //  the transition-boiling curve has its minimum at TLeidenfrost.
        //  Solves C*(exp(u)-1) - u = 0 for u = k*lDT_CHF^2 by Newton-Raphson,
        //  where C = qFilmDeriv*lDT_CHF / (2*(qCHF - qFilmLeid)).
        scalar computeTransitionDecay
        (
            const scalar qCHF,
            const scalar qFilmLeid,
            const scalar qFilmDeriv,
            const scalar lDT_CHF
        ) const;


public:

    TypeName("coolantChannelRIA");


    // Constructors

        //- Construct from patch and internal field
        coolantChannelRIAfvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        coolantChannelRIAfvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&,
            const bool valueRequired = true
        );

        //- Construct by mapping given field onto a new patch
        coolantChannelRIAfvPatchScalarField
        (
            const coolantChannelRIAfvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&,
            const bool mappingRequired = true
        );

        //- Construct as copy setting internal field reference
        coolantChannelRIAfvPatchScalarField
        (
            const coolantChannelRIAfvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new coolantChannelRIAfvPatchScalarField(*this)
            );
        }

        //- Clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new coolantChannelRIAfvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        virtual bool fixesValue() const
        {
            return true;
        }

        virtual void updateCoeffs();

        virtual void evaluate
        (
            const Pstream::commsTypes commsType = Pstream::commsTypes::blocking
        );

        virtual tmp<Field<scalar>> valueInternalCoeffs
        (
            const tmp<scalarField>&
        ) const;

        virtual tmp<Field<scalar>> valueBoundaryCoeffs
        (
            const tmp<scalarField>&
        ) const;

        virtual tmp<Field<scalar>> gradientInternalCoeffs() const;

        virtual tmp<Field<scalar>> gradientBoundaryCoeffs() const;

        virtual void write(Ostream&) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
