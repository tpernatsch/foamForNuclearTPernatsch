/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2506                                                  |
|    Copyright 2011-2016 OpenFOAM Foundation, 2017-2025 OpenCFD Ltd.          |
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
    Main object for XS parametrization and other nuclear data. Nuclear data are
    always provided in SI units and not in cm. This object generate a collection
    of ``nuclearDataOneEnergy``, one per nuclear data type (``D``,
    ``nuSigmaEff``, ...) and per energy group.

    See 2D_MSFR or 3D_SmallESFR of Godiva_SN for exhaustive examples.

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

            //- Flag to read variable XS fields. Be careful, can be intensive
            bool isReadXS_;

            //- Flag to write variable XS fields. Be careful, can be intensive
            bool isWriteXS_;

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

            //- XS variables for parametrization
            dictionary& xsVariablesDict_;

            //- Nb of XS variables
            label nxsVariables_;

            //- List variable names
            List<word> xsVariableNames_;

            //- List variable type
            List<word> xsVariableTypes_;

            //- XS variable fields, use UPtrList to avoid double free
            UPtrList<volScalarField> xsVariableFields_;

            //- XS variable fields transformed by the law provided in
            //  xsVariablesDict_
            PtrList<volScalarField> xsVariableTransformedFields_;

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

            //- Mesh displacement
            volVectorField disp_;


        // Feedbacks Coefficients

            //- Constant inverse velocity data (zoneI -> energyI -> 1/v)
            PtrList<scalarField> IVList_;

            //- Constant prompt fission emission spectrum data (zoneI -> energyI -> chi_p)
            List<PtrList<nuclearDataOneEnergy>> chiPromptList_;

            //- Constant delayed neutron emission spectrum data (zoneI -> energyI -> chi_d)
            List<PtrList<nuclearDataOneEnergy>> chiDelayedList_;

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
            label energyJ
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
        void correct();

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
