/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2406                                                  |
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
    Foam::pointKineticNeutronics

Group
    none

Description
    Run-time selectable point-kinetic based neutronics solver. It includes
    options for a time dependent external reactivity, and an external
    reactivity input coming from an FMU.

    The module solves the following point-kinetics equations:

    \f[
        \frac{dP}{dt} = \frac{\rho_{tot} - \zeta - \beta_{tot}}{\Lambda} P
            + \sum_{i=1}^N \lambda_i C_i + \frac{S}{\Lambda} \\
        \frac{dC_i}{dt} = \frac{\beta_i}{\Lambda} P
            - \lambda_i C_i \quad \forall i \in [1; N]
    \f]

    \f$ \rho_{tot} \f$: the total reactivity that can be decomposed as:
    \f[
        \rho_{tot} =
            \begin{cases}
                \text{fast}    & -\alpha_D \ln\left( \frac{T_{fuel}}{T_{fuel}^{ref}} \right) \\
                \text{thermal} &  \alpha_D \left( \sqrt{T_{fuel}} - \sqrt{T_{fuel}^{ref}} \right)
            \end{cases}
            + \alpha_{fuel} (T_{fuel} - T_{fuel}^{ref})
            + \alpha_{clad} (T_{clad} - T_{clad}^{ref})
            + \alpha_{cool}^{T} (T_{cool} - T_{cool}^{ref})
            + \alpha_{cool}^{\rho} (\rho_{cool} - \rho_{cool}^{ref}) \\
            + \alpha_{struct} (T_{struct} - T_{struct}^{ref})
            + \alpha_{struct}^{mech} (T_{struct}^{mech} - T_{struct}^{mech,ref})
            + \rho_{driveline}
            + \rho_{GEM}
            + \rho_{boron}
            + \rho_{ext}
            + \rho_{ext}^{FMU}
    \f]

    where:
    \vartable
      P           | Core power
      C_i         | Precursor power
      S           | Optional external source power
      \beta_{tot} | Total delayed neutron fraction
      \beta_i     | Delayed neutron fraction for the ith delayed neutron precursor group
      \lambda_i   | Decay constant for the ith delayed neutron precursor group
      \Lambda     | Prompt generation lifetime
      \zeta       | Optional subcritical index
      S           | Optional external neutron source power
      \alpha_D             | Doppler feedback coefficient
      \alpha_{fuel}        | Fuel axial expansion feedback coefficient
      \alpha_{clad}        | Cladding temperature feedback coefficient
      \alpha_{cool}^{T}    | Coolant temperature feedback coefficient
      \alpha_{cool}^{\rho} | Coolant density coefficient
      \alpha_{struct}      | Core radial expansion feedback coefficient
      \alpha_{struct}^{mech} | Core expansion feedback coefficient from the mechanical mesh
      \rho_{driveline}     | Control rod driveline expansion reactivity contribution
      \rho_{boron}         | Boron concentration reactivity contribution
      \rho_{ext}           | External reactivity contribution
      \rho_{ext}^{FMU}     | External reactivity contribution from an FMU
    \endtable

\vartable
    P           | Core power
    C_i         | Precursor power
    S           | Optional external source power
    \beta_{tot} | Total delayed neutron fraction
    \beta_i     | Delayed neutron fraction for the ith delayed neutron precursor group
    \lambda_i   | Decay constant for the ith delayed neutron precursor group
    \Lambda     | Prompt generation lifetime
    \zeta       | Optional subcritical index
    S           | Optional external neutron source power
    \alpha_D             | Doppler feedback coefficient
    \alpha_{fuel}        | Fuel axial expansion feedback coefficient
    \alpha_{clad}        | Cladding temperature feedback coefficient
    \alpha_{cool}^{T}    | Coolant temperature feedback coefficient
    \alpha_{cool}^{\rho} | Coolant density coefficient
    \alpha_{struct}      | Core radial expansion feedback coefficient
    \alpha_{struct}^{mech} | Core expansion feedback coefficient from the mechanical mesh
    \rho_{driveline}     | Control rod driveline expansion reactivity contribution
    \rho_{boron}         | Boron concentration reactivity contribution
    \rho_{ext}           | External reactivity contribution
    \rho_{ext}^{FMU}     | External reactivity contribution from an FMU
\endtable

Example of time dependant source entries:
\verbatim
    //- Optional time dependant external reactivity time profile
    externalReactivityTimeProfile
    {
        type        table;
        startTime   100;
        table       table
        (
            (   0   0)
            (   1   0.0005767131    ) //- i.e. 57.67131 pcm
        );
    }

    //- Optional time dependant boron reactivity time profile
    boronReactivityTimeProfile
    {
        type        table;
        table       table
        (
            (   100   0   )
            (   101   1e-5   )
        );
    }
\endverbatim

Optional external reactivity input from FMU in nuclearData file.
The default initial value is 0.
\verbatim
    externalReactivityNameFromFMU   gfExtReact;
\endverbatim

Optional external neutron source file "externalSource"
\verbatim
    //- If true, add the external source to the point-kinetics equations
    isExternalSource    true;

    //- External source mode
    externalSourceMode  transient; // powerMonitoring // transient

    //- Beam energy per source particle
    beamEnergy          1.28160e-10; // [J/proton], 800 MeV/proton

    //- Neutron yield of the spallation of proton on heavy metal
    nuSource            15.5; // [neutron/proton]

    //- Optional external source input from FMU in externalSource file.
    //- The externalSourceMode must be set to transient. The default initial value
    //-is 1.
    externalSourceModulationNameFromFMU     gfExternalSourceMod;

    //- External modulation time profile (time, modulation factor)
    externalSourceModulationTimeProfile
    {
        type        table;
        startTime   150;
        table       table
        (
            (0      1)
            (20     2)
        );
    }
\endverbatim

\mainauthor of this file (and associated .C or included .H files):
    Stefan Radman <stefanradman92@gmail.com; stefan.radman@epfl.ch;>;
    Carlo Fiorina <carlo.fiorina@outlook.com; carlo.fiorina@epfl.ch;>, EPFL (Switzerland);
    Thomas Guilbaud <thomas.guilbaud@epfl.ch>, EPFL/Transmutex SA (Switzerland)

SourceFiles
    pointKineticNeutronics.C

\*---------------------------------------------------------------------------*/

#ifndef pointKineticNeutronics_H
#define pointKineticNeutronics_H

#include "neutronics.H"
#include "Function1.H"
#include "interpolationTable.H"
#include "timeProfile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers{

/*---------------------------------------------------------------------------*\
                          Class pointKineticNeutronics Declaration
\*---------------------------------------------------------------------------*/

class pointKineticNeutronics
:
    public neutronics
{
protected:

    //-
    IOdictionary nuclearData_;

    //- Total power, sum of fission power and decay power
    scalar power_;

    //- Point kinetics operates on fissionPower
    scalar fissionPower_;

    //- Decay power value set by decayPowerPtr_
    scalar decayPower_;

    //- Optional decay power
    // autoPtr<Function1<scalar>> decayPowerPtr_;
    timeProfile decayPowerTimeProfile_;

    //- Time offset for the calcuation of the decayPower based on the
    //  provided Function1, if any
    // scalar decayPowerStartTime_;

    //- Due to whatever the user might please
    timeProfile externalReactivityTimeProfile_;

    //- External source mode
    //  - transient: user input power variation by modulation
    //  - powerMonitoring: adapt source power to keep subcritical power constant
    string externalSourceMode_;

    //- Anti-reactivity due to the subcritical configuration of the core. Also
    //  Noted as zeta
    scalar subcriticalIndex_;

    //- ksrc of the system. Allows to compute the subcritical-index in case it
    //  is not provided
    scalar ksrc_;

    //- Avg number of neutrons generated by fission
    scalar nuFission_;
    //- Avg number of neutrons generated by high energy reaction (e.g. spallation)
    scalar nuSource_;
    //- Avg energy released per fission
    scalar energyPerFission_;
    //- Energy of the particle beam in case external source from particle
    //  accelerator
    scalar beamEnergy_;

    // G_0 = energyPerFission*nuSource/(nuFission*beamEnergy)
    scalar intrinsicGain_;

    //- Source modulation value set by externalSourceModulationTimeProfile_ (m_I)
    scalar externalSourceModulation_;

    //- External source power (S)
    scalar externalSourcePower_;
    scalar externalSourcePowerOld_;
    scalar externalSourcePowerRef_;

    //- External source beam power (P_src)
    scalar externalSourceBeamPower_;

    //- External source beam intensity (I_src)
    scalar externalSourceBeamIntensity_;

    //- Power to target during subcritical powerMonitoring mode using the
    //  modulation of the external source
    scalar powerTarget_;

    //- Record of the subcritical power in case of the powerMonitoring mode.
    List<Tuple2<scalar, scalar>> powerRecordPtr_;

    //- Modulation factor to change the efficienty of the source power
    //  modulation (k_m)
    scalar modulationFactor_;

    //- Delay between the power of the subcritical system and the external
    //  source modulation
    scalar powerModulationDelay_;

    //- Used to set an initial steady state for liquid fuel
    scalar precEquilibriumReactivity_;

    //- Estimated value of Beta for liquid fuel
    scalar liquidFuelBeta_;

    //-
    scalar totalReactivity_;

    //-
    scalar promptGenerationTime_;

    //- Total effective delayed neutron fraction
    scalar beta_;

    //- Effective delayed fractions for each delayed group
    scalarList betas_;

    //- Decay constants for each delayed group
    scalarList lambdas_;

    //- Precursor "power" for each delayed group in W to be compatible with
    //  point kinetics being solved for power rather than neutron density
    scalarList precursorPowers_;

    //- Number of delayed groups
    label delayedGroups_;

    //- Old time power, precursorPowers to be used to calculate derivatives,
    //  as well as timeIndex to check if time step has changed (can't think
    //  of a better way to get something equivalent to the finalInnerIter flag)
    label timeIndex_;
    scalar powerOld_;
    scalar fissionPowerOld_;
    scalarList precursorPowersOld_;

    // default precursor field. Only used for liquid fuel
    volScalarField defaultPrec_;

    //- Feedback coefficients in 1/K (i.e. 1e5 pcm/K)
    scalar coeffTFuel_;
    scalar coeffTClad_;
    scalar coeffTCool_;
    scalar coeffRhoCool_;
    scalar coeffTStruct_;
    scalar coeffTStructMech_;

    //- Absolute driveline expansion coefficient in m/K
    scalar coeffDrivelineExp_;

    //- Optional Boron reactivity curve
    timeProfile boronReactivityTimeProfile_;

    //- Reference integral values for feedbacks
    scalar TFuelRef_;
    scalar TCladRef_;
    scalar TCoolRef_;
    scalar rhoCoolRef_;
    scalar TStructRef_;
    scalar TStructMechRef_;
    scalar TDrivelineRef_;

    //- Ptrs to coupling fields (i.e. defined on their original mesh)
    const volScalarField* TFuelOrig_;
    const volScalarField* TCladOrig_;
    const volScalarField* TCoolOrig_;
    const volScalarField* rhoCoolOrig_;
    const volScalarField* TStructOrig_;
    const volScalarField* TStructMechOrig_;
    const volScalarField* powerDensityOrig_;
    const volScalarField* powerDensityToLiquidOrig_;
    const volVectorField* UOrig_;
    const volScalarField* alphaOrig_;
    const volScalarField* alphatOrig_;
    const volScalarField* muOrig_;

    //- Mapped fields from TH
    volScalarField TFuel_;
    volScalarField TClad_;
    volScalarField TCool_;
    volScalarField rhoCool_;
    volScalarField TStruct_; //- This is also used to compute TDriveline
    //- Mapped fields from TM
    volScalarField TStructMech_;

    //- Also mapped, but these are autoPtr as these field are not needed at all
    //  unless doing liquidFuel. Thus, avoid wasting memory by only setting
    //  these fields if they are actually needed
    autoPtr<volVectorField> UPtr_;
    autoPtr<volScalarField> alphaPtr_;
    autoPtr<volScalarField> alphatPtr_;
    autoPtr<volScalarField> muPtr_;

    //- Fields derived from mapped fields
    autoPtr<surfaceScalarField> phiPtr_;
    autoPtr<volScalarField> diffCoeffPrecPtr_;

    //- Fields necessary when an albedo boundary condition has been used
    //  to calculate the initial fluxes
    volScalarField Dalbedo_;
    volScalarField fluxStarAlbedo_;

    //-
    volScalarField oneGroupFlux_;

    //- Initial flux
    volScalarField initOneGroupFlux_;

    //- Normalized initial flux
    volScalarField initOneGroupFluxN_;

    //- Integral of the square of the normalized flux
    scalar domainIntegratedInitOneGroupFluxN_;

    //- Number of energy groups (for scaling fluxes of a possibly existing
    //  previous simulation)
    label energyGroups_;

    //- Group-wise fluxes
    PtrList<volScalarField> fluxes_;

    //- Delayed-group-wise precursor concentrations
    PtrList<volScalarField> precursors_;

    //- Delayed-group-wise precursor concentrations resulting from point
    //  kinetics calculations with liquid fuel. This second set of fields is
    //  necessary because the unit for precursors in point kinetics is power,
    //   instead of 1/m3
    PtrList<volScalarField> precPK_;//reconstructed based on porosity
    PtrList<volScalarField> precPKStar_;//solved for

    //- If true, fuel temperature feedback (Doppler) computed as
    //  deltaReactivity = -coeffFastDoppler*ln(TFuel/TFuelRef) with
    //  coeffFastDoppler being the Doppler constant. If false, Doppler feedback
    //  computed as deltaReactivity = coeffTFuel*(TFuel-TFuelRef)
    bool fastNeutrons_;

    //- Doppler constant if doing fastNeutrons
    scalar coeffFastDoppler_;

    //- Fields to mark which cells are part of X and which aren't, used for
    //  computing the 0-D feedback values (TFuel, TStructures, etc...).
    //  Indeed, the spatial extent of these marker fields should reflect
    //  the spatial extent of the domain that was perturbed in neutronics
    //  calculations to compute feedback coefficients

    //- For TFuel and TCladding
    volScalarField fuelFeedbackCellField_;

    //- For rhoCoolant and TCoolant
    volScalarField coolFeedbackCellField_;

    //- For TStructures
    volScalarField structFeedbackCellField_;

    //- For TStructures from TM
    volScalarField structMechFeedbackCellField_;

    //- For the control rod driveline
    volScalarField drivelineFeedbackCellField_;

    //- Pairs of control rod insertion value (in m w.r.t. a reference) vs
    // reactivity in absolute units (i.e. 1e5 pcm)
    List<Pair<scalar>> controlRodReactivityMap_;

    //- Schmidt number. Only used for liquid fuel
    scalar ScNo_;

    //- Flag to initialize the concentration of precursors to their
    //  equilibtium value at the beginning of the simulation
    bool initPrecursorsLiquidFuel_;

    //- GEM related

    //- Reactivity map
    List<Pair<scalar>> GEMReactivityMap_;

    //- Reference height of sodium in GEM
    scalar GEMSodiumLevelRef_;

    //- Ref to volumetric flux on fluid mesh
    const surfaceScalarField* phiOrig_;

    //- Reference surface-integrated value of phi
    scalar intPhiRef_;

    //- List of face labels over which to compute integral flow values
    labelList phiFaces_;

    //- Surface areas of the phiFaces_ (I cannot access fvMesh from polyMesh,
    //  which is the only mesh-related object I can access from the mapping,
    //  thus magSf needs to be computed by calling the faceArea() method
    //  of polyMesh. To avoid that, I calculated the areas only once and
    //  store them in pgiMagSf_)
    scalarList phiMagSf_;

    scalar pTotOld_;

    scalar residual_;

private:

    //- Disallow default bitwise copy construct
    pointKineticNeutronics(const pointKineticNeutronics&);

    //- Disallow default bitwise assignment
    void operator=(const pointKineticNeutronics&);

public:

    //- Runtime type information
        TypeName("pointKinetics");

    //- Constructors

        //- Construct from mesh, read data from IOdictionaries
        pointKineticNeutronics(fvMesh& mesh);

    //- Destructor
        virtual ~pointKineticNeutronics(){}

    //- Member Functions

        //- Calc keff from reactivity
        virtual scalar keff() const override
        {
            return 1.0/(1.0-totalReactivity_);
        }

        //- Update all the oneGroupFlux related quantities following changes in
        //  the oneGroupFlux after proper initialization
        void setInitOneGroupFlux();

        //- Used to set the *FeedbackCellField fields
        void setFeedbackCellField
        (
            volScalarField& feedbackCellField,
            const word& keyword
        );

        //- Set beam parameter from the external source power S
        void setBeamParameters();

        //- Compute driveline expansion reactivity from controlRod reactivity
        //  map given the absolute driveline expansion in meters
        scalar calcDrivelineReactivity(const scalar& drivelineExpansion);

        //- Compute GEM sodium height and their reactivity contributions
        Pair<scalar> calcGEMLevelAndReactivity();

        //- Solve pointKinetics
        virtual void correctPhysics();
        virtual void correctTightlyCoupledPhysics();

        virtual scalar maxDeltaT();

        virtual scalar getResidual()
        {
            return residual_;
        }

        //- To split strings. The file where I got this function originally
        //  defined is in classes/thermalHydraulics/models/include. However,
        //  including things from thermalHydraulics in neutronics, even if it's
        //  just a helper function such as this one, is tremendously
        //  un-elegant. I don't want the code structure to be spaghetti, so
        //  I am re-implementing the function here

        //- The deliminter is a char, so it must be provided in single-
        //  quotation marks. Type can be either word or string
        template<class Type>
        List<Type> split(const Type& compoundItem, char d)
        {
            List<Type> items(1, "");
            int i = 0;
            for (char const &c : compoundItem)
            {
                if (c == d)
                {
                    i += 1;
                    items.append("");
                }
                else items[i] += c;
            }
            return items;
        }

        //- Access

            //- Return fluxes
            virtual const PtrList<volScalarField> fluxes() const
            {
                return fluxes_;
            };

            //- Return precursors
            virtual const PtrList<volScalarField> precursors() const
            {
                return precursors_;
            };

            //- Return one group flux
            virtual const volScalarField& oneGroupFlux() const
            {
                return oneGroupFlux_;
            }

            //- Return power
            virtual scalar power() const
            {
                return power_;
            };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
