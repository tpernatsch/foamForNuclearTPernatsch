/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2412                                                  |
|    Copyright 2011-2016 OpenFOAM Foundation, 2017-2024 OpenCFD Ltd.          |
-------------------------------------------------------------------------------
License
    This file is part of GeN-Foam.

    GeN-Foam is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    GeN-Foam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    This offering is not approved or endorsed by the OpenFOAM Foundation nor
    OpenCFD Limited, producer and distributor of the OpenFOAM(R)software via
    www.openfoam.com, and owner of the OPENFOAM(R) and OpenCFD(R) trademarks.

    This particular snippet of code is developed according to the developer's
    knowledge and experience in OpenFOAM. The users should be aware that
    there is a chance of bugs in the code, though we've thoroughly test it.
    The source code may not be in the OpenFOAM coding style, and it might not
    be making use of inheritance of classes to full extent.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::XS

Description
    Class for XS parametrization and other nuclear data. Nuclear data are always
    provided in MKSA units and not in cm.

See 2D_MSFR or 3D_SmallESFR of Godiva_SN for exhaustive examples.

<div class="groupHeader">
\table
    Parameter           | Type      | Req'd | Default | Description
    energyGroups        | Integer   | No    | 1       | Number of energy groups
    precGroups          | Integer   | No    | 1       | Number of precursor groups
    legendreMoments     | Integer   | No    | 0       | Number of Legendre moments
    fastNeutrons        | Boolean   | No    | True    | Flag for fast neutron spectrum, if true, use log(Tfuel) law, else use sqrt(Tfuel) law
    adjustDiscFactors   | Boolean   | No    | False   | Flag to apply the discontinuity factor adjustement
    useGivenDiscFactors | Boolean   | No    | False   | Flag to use homogeneous discontinuity factors provided in state/zone
    groupsWoDF          | List<int> | No    | Empty   | List of energy group index to not apply discontinuity factors
    doNotParametrize    | List<int> | No    | Empty   | List of energy group index to not use parametrization (use only the values from the reference state)
    ScNo                | Scalar    | No    | 1       | Schmidt number for diffusion of precursors
    polyharmonicSplineMode | Label  | No    | 1 (linear) | Polyharmonic spline function mode (1: \f$ r \f$; 2: \f$ r^2 ln(r) \f$; 3: \f$ r^3 \f$; 4: \f$ r^4 ln(r) \f$)
    axialOrientation    | Vector    | No    | (0,0,1) | Fuel axial orientation
    isLowMemory         | Boolean   | No    | False   | Low memory flag, useful for SN solver
\endtable
</div>

The following table lists the possible perturb state keywords to be used in the
state dictionaries. it is possible to use multiple parameters in the same state
(e.g Tfuel and rhoCool).

<div class="groupHeader">
\table
    Parameter   | Type   | Req'd | Default | Description
    Tfuel       | Scalar | No    | 0.01    | Fuel temperature in K
    Tclad       | Scalar | No    | 0       | Cladding temperature in K
    Tcool       | Scalar | No    | 0       | Coolant temperature in K
    TstructMech | Scalar | No    | 0       | Mechanical structure temperature in K
    rhoCool     | Scalar | No    | 0       | Coolant density in kg/m3
    axExp       | Scalar | No    | 0       | Axial expansion from nominal
    radExp      | Scalar | No    | 0       | Radial expansion from nominal
\endtable
</div>

Here follow the list of group wise nuclear data for each cellzone to be provided
for GeN-Foam. Each subscript \f$g\f$ refers the gth energy group. Each \f$k\f$
refers to the kth precursor group.

<div class="groupHeader">
\table
    GeN-Foam Parameter  | Param    | Req'd | Default | Notation                 | Description
    fuelFraction        | Constant | No    | 1       | \f$ \alpha_{fuel} \f$    | Volume of fuel per lattice volume
    secondaryPowerVolumeFraction | Constant | No | 1 |                          | Volume fraction of secondary power-producing structure e.g., graphite in MSRs
    fractionToSecondaryPower | Constant | No | 0     |                          | Fraction of total power that goes to secondary power-producing structure
    IV                  | Constant | Yes   | -       | \f$ 1/v_g \f$            | Inverse velocity
    D                   | Variable | Yes   | -       | \f$ D_g \f$              | Diffusion coefficient
    nuSigmaEff          | Variable | Yes   | -       | \f$ \nu \Sigma_{f,eff, g} \f$ | Fission neutrons cross-section
    sigmaPow            | Variable | Yes   | -       | \f$ \kappa \Sigma_{f, g} \f$ | Fission energy release cross-section
    sigmaRemoval        | Variable | Yes   | -       | \f$ \Sigma_{r, g} \f$ | Removal cross-section (\f$\Sigma_{abs} + \Sigma_{capture} + \Sigma_{s, g \rightarrow g}\f$)
    scatteringMatrixP0  | Variable | Yes   | -       | \f$ \Sigma_{s, g \rightarrow g'} \f$ | Scattering cross-section matrix
    chiPrompt           | Constant | Yes   | -       | \f$ \chi_{p, g} \f$      | Prompt fission emission spectrum
    chiDelayed          | Constant | Yes   | -       | \f$ \chi_{d, g} \f$      | Delayed fission spectrum
    Beta                | Constant | Yes   | -       | \f$ \beta_k \f$          | Delayed neutron fraction
    BetaTot             | Constant | Yes   | -       | \f$ \beta_{tot} \f$      | Total delayed neutron fraction
    lambda              | Constant | Yes   | -       | \f$ \lambda_k \f$        | Precursor decay constant
    discFactor          | Constant | Yes   | -       | \f$ \gamma_g \f$         | Discontinuity factors
    integralFlux        | Constant | Yes   | -       | \f$ \Phi_g \f$           | Integral flux for adapting disc factors
\endtable
</div>

Example of dictionary entries for perturbed state XS in the nuclearData file:
\verbatim
states
(
    reference // Mandatory state
    {
        // Reference parameters
        Tfuel       1000;
        rhoCool     1000;
        // Assume the other parameters are not used

        zones
        (
            core
            {
                fuelFraction    ...;
                IV  ...;
                ...
            }
            reflector
            {
                ...
            }
        );
    }
    Tfuel1200K // Additional state, arbitrary state name
    {
        // Perturb only Tfuel
        Tfuel       1200;
        // Assume by default that rhoCool is not perturb and uses rhoCool = 1000

        zones
        (
            core
            {
                ...
            }
        );
    }
    TfuelAndRhoHot
    {
        // Use both Tfuel and rhoCool parameters to perturb the state
        Tfuel       1200;
        rhoCool     900;

        zones
        (
            core
            {
                ...
            }
            reflector
            {
                ...
            }
        );
    }
);
\endverbatim

\mainauthor of this file (and associated .C or included .H files)
    Carlo Fiorina <carlo.fiorina@outlook.com; carlo.fiorina@epfl.ch;>, EPFL (Switzerland);
    Thomas Guilbaud <thomas.guilbaud@epfl.ch>, EPFL/Transmutex SA (Switzerland)

SourceFiles
    XS.C

\*---------------------------------------------------------------------------*/

#ifndef XS_H
#define XS_H

#include "volFields.H"
#include "rhoThermo.H"
#include "fvMatrices.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "coordinateSystem.H"
#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "volPointInterpolation.H"
#include "meshToMesh.H"
#include "nuclearDataOneEnergy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                         Class thermoMechanicsSolver Declaration
\*---------------------------------------------------------------------------*/

class XS
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        XS(const XS&);

        //- Disallow default bitwise assignment
        void operator=(const XS&);

protected:

    // Protected data

        //- Reference to mesh
            const fvMesh& mesh_;

        // Data

            //- Nuclear data dictionary
            IOdictionary nuclearData_;

            //- Flag for low memory consumption when creating XS field
            bool isLowMemory_;

            //- Number of energy groups
            label energyGroups_;

            //- Number of precursor groups
            label precGroups_;

            //- Number of Legendre moments
            label legendreMoments_;

            //- Fuel axial orientation vector, default along the Z-axis
            vector axialOrientation_;

            //- Schmidt number for diffusion of precursors
            scalar ScNo_;

            //- Polyharmonic spline function mode
            label polyharmonicSplineMode_;

            //- List of perturbed reactor states
            PtrList<entry> states_;

            //- Reference to the reference reactor state
            dictionary& referenceState_;

            //- List of all cellzone in the reference reactor state
            PtrList<entry> referenceZones_;

            //- Total number of cellzone
            label zoneNumber_;

        // Data fields

            //- Energy group-wise list of inverse velocity field
            PtrList<volScalarField> IV_;

            //- Energy group-wise list of diffusion coefficient field
            PtrList<volScalarField> D_;

            //- Energy group-wise list of fission neutrons cross-section field
            PtrList<volScalarField> nuSigmaEff_;

            //- Energy group-wise list of fission energy released cross-section field
            PtrList<volScalarField> sigmaPow_;

            //- Energy group-wise list of disapearance cross-section field
            PtrList<volScalarField> sigmaRemoval_;

            //- Moments/EnergyJ/EnergyI array of scattering matrices cross-section field
            List<List<PtrList<volScalarField>>> sigmaFromTo_;

            //- Energy group-wise list of prompt fission emission spectrum field
            PtrList<volScalarField> chiPrompt_;

            //- Energy group-wise list of delayed fission spectrum field
            PtrList<volScalarField> chiDelayed_;

            //- Precursor group-wise list of delayed neutron fraction field
            PtrList<volScalarField> Beta_;

            //- Total delayed neutron fraction field
            volScalarField BetaTot_;

            //- Precursor group-wise list of delayed neutron precursor decay constant field
            PtrList<volScalarField> lambda_;

            //- Fuel fraction field
            volScalarField fuelFraction_;

            //- Volume fraction of secondary power-producing structure e.g., graphite in MSRs
            volScalarField secondaryPowerVolumeFraction_;

            //- Fraction of total power that goes to secondary power-producing structure
            volScalarField fractionToSecondaryPower_;

            //- Energy group-wise list of discontinuity factor field
            PtrList<volScalarField> discFactor_;

            //- 1 if there is transfer from a group I to J in any cellZone if the XS is higher than 5.0e-6
            SquareMatrix<label> sigmaFromToYesNo_;


        // Disc factors

            //- Fuel fraction used to determine vol power IN the fuel starting from vol power from sigmaPow
            PtrList<scalar> fuelFractionList_;

            //- Volume fraction of secondary power-producing structure e.g., graphite in MSRs
            PtrList<scalar> secondaryPowerVolumeFractionList_;

            //- Fraction of total power that goes to secondary power-producing structure
            PtrList<scalar> fractionToSecondaryPowerList_;

            //- Adjust or not discontinuity factors
            PtrList<bool> dfAdjustList_;

            //- Energy group-wise discontinuity factors field
            PtrList<scalarField> discFactorList_;

            //- Energy group-wise integral fluxes field for adapting discontinuity factors
            PtrList<scalarField> integralFluxList_;

            //- Flag to ajust the discontinuity factors
            bool adjustDiscFactors_;

            //- Flag to use homogeneous discontinuity factors provided in state/zone from dictionary
            bool useGivenDiscFactors_;

            //- List of energy group index to not apply discontinuity factors
            List<int> groupsWoDF_;


        // Parametrization flags

            //- Fast spectrum flag. If true, the cross-section dependence w.r.t
            //- the fuel temperature follows a log(T) law. Otherwise follows a
            //- sqrt(T) law for thermal spectrum.
            bool fastNeutrons_;

            //- List of energy group index to not use parametrization (use only the values from the reference state)
            List<int> doNotParametrize_;


        // Fields

            //- Radial expansion field
            volScalarField radExp_;

            //- Axial expansion field
            volScalarField axExp_;

            //- Logarithmic fuel temperature field
            volScalarField logT_;

            //- Square root fuel temperature field
            volScalarField sqrtT_;


        // Feedbacks Coefficients

            //- Constant inverse velocity data (zoneI -> energyI -> 1/v)
            PtrList<scalarField> IVList_;

            //- Constant prompt fission emission spectrum data (zoneI -> energyI -> chi_p)
            PtrList<scalarField> chiPromptList_;

            //- Constant delayed neutron emission spectrum data (zoneI -> energyI -> chi_d)
            PtrList<scalarField> chiDelayedList_;

            //- Constant delayed neutron fraction data (zoneI -> precursorK -> beta)
            PtrList<scalarField> BetaList_;

            //- Constant total delayed neutron fraction data (zoneI -> beta_tot)
            PtrList<scalar> BetaTotList_;

            //- Constant precursor decay constant data (zoneI -> precursorK -> lambda)
            PtrList<scalarField> lambdaList_;

            //- Variable diffusion coefficient data (zoneI -> energyI -> XS)
            List<PtrList<nuclearDataOneEnergy>> DList_;

            //- Variable fission neutron XS data (zoneI -> energyI -> XS)
            List<PtrList<nuclearDataOneEnergy>> nuSigmaEffList_;

            //- Variable fission energy release XS data (zoneI -> energyI -> XS)
            List<PtrList<nuclearDataOneEnergy>> sigmaPowList_;

            //- Variable removal XS data (zoneI -> energyI -> XS)
            List<PtrList<nuclearDataOneEnergy>> sigmaRemovalList_;

            //- Variable scateering matrix XS data (zoneI -> momentI -> energyJ -> energyI -> XS)
            List<List<List<PtrList<nuclearDataOneEnergy>>>> sigmaFromToList_;


        // For control rods

            //- Control rod dictionary
            IOdictionary CRmove_;
            PtrList<entry> CRentries_;
            label CRNumber_;
            PtrList<scalar> CRstart_;
            PtrList<scalar> CRfinish_;
            PtrList<scalar> CRspeed_;
            PtrList<word>   CRFollowerName_;
            PtrList<scalar> CRinitialPosition_;
            PtrList<scalar> CRposition_;
            PtrList<scalar> initialDistanceFromMeshCR_;

public:

    //- Runtime type information
    TypeName("XS");

    // Declare run-time constructor selection table

    declareRunTimeSelectionTable
    (
        autoPtr,
        XS,
        dictionary,
        (const fvMesh& mesh),
        (mesh)
    );

    // Access

        label energyGroups() const
        {
            return energyGroups_;
        }
        label precGroups() const
        {
            return precGroups_;
        }
        label legendreMoments() const
        {
            return legendreMoments_;
        }
        const PtrList<volScalarField> IV() const
        {
            return IV_;
        }
        const PtrList<volScalarField> D() const
        {
            return D_;
        }
        const PtrList<volScalarField> nuSigmaEff() const
        {
            return nuSigmaEff_;
        }
        const PtrList<volScalarField> sigmaPow() const
        {
            return sigmaPow_;
        }
        const PtrList<volScalarField> sigmaRemoval() const
        {
            return sigmaRemoval_;
        }
        const List<List<PtrList<volScalarField>>> sigmaFromTo() const
        {
            return sigmaFromTo_;
        }

        //- Create sigmaFromTo field on the fly to avoid storing all of them in
        //- memory (N_moment * N_energy * N_energy)
        Foam::tmp<Foam::volScalarField> sigmaFromTo
        (
            label momentI,
            label energyI,
            label energyJ,
            const volScalarField& Tfuel,
            const volScalarField& Tclad,
            const volScalarField& rhoCool,
            const volScalarField& Tcool,
            const volVectorField& Disp,
            const volScalarField& TStructMech
        );

        const PtrList<volScalarField> chiPrompt() const
        {
            return chiPrompt_;
        }
        const PtrList<volScalarField> chiDelayed() const
        {
            return chiDelayed_;
        }
        const PtrList<volScalarField> Beta() const
        {
            return Beta_;
        }
        const volScalarField BetaTot() const
        {
            return BetaTot_;
        }
        const PtrList<volScalarField> lambda() const
        {
            return lambda_;
        }
        const volScalarField fuelFraction() const
        {
            return fuelFraction_;
        }
        const volScalarField secondaryPowerVolumeFraction() const
        {
            return secondaryPowerVolumeFraction_;
        }
        const volScalarField fractionToSecondaryPower() const
        {
            return fractionToSecondaryPower_;
        }
        const PtrList<volScalarField> discFactor() const
        {
            return discFactor_;
        }
        scalar ScNo() const
        {
            return ScNo_;
        }
        const SquareMatrix<label>& sigmaFromToYesNo() const
        {
            return sigmaFromToYesNo_;
        }

    // Constructors

        //- Construct from mesh and thermo-mechanical properties
        XS
        (
            const fvMesh& mesh
        );

        //- Return a pointer to a new XS
        static autoPtr<XS> New
        (
            const fvMesh& mesh
        );

    //- Destructor
        virtual ~XS();

    // Member Functions

        //- Update XS
        void correct
        (
            const volScalarField& Tfuel,
            const volScalarField& Tclad,
            const volScalarField& rhoCool,
            const volScalarField& Tcool,
            const volVectorField& Disp,
            const volScalarField& TStructMech
        );

        //- Initialize the XS data set
        void init();

        //-
        void adjustDiscFactors(const PtrList<volScalarField>& fluxStar);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
