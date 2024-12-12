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
    Foam::diffusionNeutronics

Description
    Subclass for diffusion neutronics.

The module solves the multi-group diffusion equations:

\f[
    \frac{1}{v_g} \frac{\partial \phi_g}{\partial t} =
        \nabla D_g \nabla \phi_g
        + \frac{\chi_{p,g}}{k_{eff}} (1-\beta_{tot}) \nu \Sigma_{f,g} \phi_g
        - \Sigma_{r,g} \phi_g
        + \frac{\chi_{p,g}}{k_{eff}} (1-\beta_{tot}) S_{n,g}
        + \chi_{d,g} S_d
        + S_{s,g}
        + S_g
        \quad \forall g \in [1; G] \\
    \frac{\partial C_k}{\partial t} + \nabla \cdot (\boldsymbol{u} C_k) =
        \frac{\beta_k}{k_{eff}} \sum_h \nu \Sigma_{f,h} \phi_h
        - \lambda_k C_k
        \quad \forall k \in [1; N]
\f]

with the explicit source term:

\f[
    S_{n,g} = \sum_{h \neq g} \nu \Sigma_{f,h} \phi_h \\
    S_d = \sum_k \lambda_k C_k \\
    S_{s,g} = \sum_{h \neq g} \Sigma_{h \rightarrow g} \phi_h
\f]

where:
\vartable
    \phi_g         | Neutron flux for the gth energy group
    C_k            | Concentration of the kth delayed neutron precursor group
    D_g            | Neutron diffusion coefficient for the gth energy group
    \chi_{p,g}     | Prompt neutron yield for the gth energy group
    \chi_{d,g}     | Delayed neutron yield for the ith energy group
    k_{eff}        | Effective multiplication factor
    \beta_{tot}    | Total delayed neutron fraction
    \beta_k        | Delayed neutron fraction for the kth delayed neutron precursor group
    \nu            | Average number of neutrons per fission
    \Sigma_{f,g}   | Fission cross section for the gth energy group
    \Sigma_{r,g}   | Removal (disappearance) cross section for the gth energy group
    \Sigma_{h \rightarrow g} | Group-transfer cross section from the hth to the gth energy group
    S_{n,g}        | Fission neutron source from neutron energy groups others than the gth
    S_d            | Delayed neutron source
    S_{s,g}        | Scattering neutron source from neutron energy groups others than the gth
    S_g            | External neutron source for the gth energy group
    \boldsymbol{u} | Fuel velocity
    v_g            | Average neutron velocity for the gth energy group
    \lambda_k      | Decay constant for the kth delayed neutron precursor group
\endvartable


In the neutronicsProperties file, for keff calculation in a critical
reactor:
\verbatim
    model                       diffusionNeutronics;
    eigenvalueNeutronics        true;
\endverbatim

For an subcritical reactor driven by an external neutron source:
\verbatim
    model                       diffusionNeutronics;
    eigenvalueNeutronics        false;
    externalSourceNeutronics    true;
\endverbatim

<b>Reference publications</b>

<ol type="1">
    <li> Carlo Fiorina, Ivor Clifford, Manuele Aufiero, Konstantin Mikityuk,
    GeN-Foam: a novel OpenFOAM® based multi-physics solver for 2D/3D transient
    analysis of nuclear reactors, Nuclear Engineering and Design, Volume 294,
    2015, Pages 24-37, ISSN 0029-5493, https://doi.org/10.1016/j.nucengdes.2015.05.035.

    <li> Carlo Fiorina, Konstantin Mikityuk, "Application of the new GeN-Foam
    multi-physics solver to the European Sodium Fast Reactor and verification
    against available codes", Proceedings of ICAPP 2015, May 03-06, 2015 -
    Nice (France), Paper 15226

    <li> Carlo Fiorina, Nordine Kerkar, Konstantin Mikityuk, Pablo Rubiolo, Andreas
    Pautz, "Development and verification of the neutron diffusion solver for the
    GeN-Foam multi-physics platform", Annals of Nuclear Energy, 96, 2016,
    pages 212-222, DOI: http://dx.doi.org/10.1016/j.anucene.2016.05.023.
</ol>

\mainauthor of this file (and associated .C or included .H files):
    Carlo Fiorina <carlo.fiorina@outlook.com; carlo.fiorina@epfl.ch;>, EPFL (Switzerland);
    Thomas Guilbaud <thomas.guilbaud@epfl.ch>, EPFL/Transmutex SA (Switzerland)

SourceFiles
    diffusionNeutronics.C

\*---------------------------------------------------------------------------*/

#ifndef diffusionNeutronics_H
#define diffusionNeutronics_H

#include "neutronics.H"
#include "XS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class diffusionNeutronics Declaration
\*---------------------------------------------------------------------------*/

class diffusionNeutronics
:
    public neutronics
{
protected:

    //- Cross sections
    XS xs_;


    //- Fields
    volScalarField Dalbedo_;
    PtrList<volScalarField> flux_;
    PtrList<volScalarField> fluxStar_;
    PtrList<volScalarField> externalSourceFlux_;
    PtrList<volScalarField> prec_; // prec per total volume
    PtrList<volScalarField> precStar_; // prec per volume salt
    volScalarField fluxStarAlbedo_;
    volScalarField defaultFlux_;
    volScalarField defaultExternalSourceFlux_;
    volScalarField defaultPrec_;

    //- Derived fields
    // volScalarField oneGroupFlux_;
    volScalarField neutroSource_;
    volScalarField delayedNeutroSource_;
    volScalarField scatteringSourceExtra_;

    //- Ptrs to coupling fields (i.e. defined on their original mesh)
    // const volScalarField* TFuelOrig_;
    // const volScalarField* TCladOrig_;
    // const volScalarField* TCoolOrig_;
    // const volScalarField* rhoCoolOrig_;
    // const volScalarField* TStructMechOrig_;
    // const volVectorField* UOrig_;
    // const volScalarField* alphaOrig_; //- Total volume fraction of fluid(s)
    // const volScalarField* alphatOrig_; //- Thermal diffusivity in kg/m/s
    // const volScalarField* muOrig_; //- Molecular viscosity in kg/m/s (== Pa*s)

    //- Mapped fields from TH
    // volScalarField TFuel_;
    // volScalarField TClad_;
    // volScalarField TCool_;
    // volScalarField rhoCool_;
    // //- Mapped fields from TM
    // volScalarField TStructMech_;

    //- Also mapped, but these are autoPtr as these field are not needed at all
    //  unless doing liquidFuel. Thus, avoid wasting memory by only setting
    //  these fields if they are actually needed
    // autoPtr<volVectorField> UPtr_;
    // autoPtr<volScalarField> alphaPtr_;
    // autoPtr<volScalarField> alphatPtr_;
    // autoPtr<volScalarField> muPtr_;

    //- Fields derived from mapped fields
    // autoPtr<surfaceScalarField> phiPtr_;
    // autoPtr<volScalarField> diffCoeffPrecPtr_;

    // scalar pTotOld_;

    // scalar residual_;

    //- multiRegion
    // bool multiRegion_;


private:

    //- Disallow default bitwise copy construct
    diffusionNeutronics(const diffusionNeutronics&);

    //- Disallow default bitwise assignment
    void operator=(const diffusionNeutronics&);

public:

    //- Runtime type information
        TypeName("diffusionNeutronics");

    //- Constructors
        diffusionNeutronics(fvMesh& mesh);

    //- Destructor
        virtual ~diffusionNeutronics();

    //- Member Functions

        //- Return fluxes
        virtual const PtrList<volScalarField> fluxes() const
        {
            return flux_;
        }

        //- Return precursors
        virtual const PtrList<volScalarField> precursors() const
        {
            return prec_;
        }

        //- Return one group flux
        // virtual const volScalarField& oneGroupFlux() const
        // {
        //     return oneGroupFlux_;
        // }

        //- Return power
        virtual scalar power() const
        {
            scalar pTot = 0;

            forAll(flux_, energyI)
            {
                pTot += fvc::domainIntegrate
                (
                    flux_[energyI] * xs_.sigmaPow()[energyI]
                ).value();
            }
            return pTot;
        }

        //- Correct/update the properties
        virtual void correctPhysics();
        virtual void correctTightlyCoupledPhysics();

        virtual scalar maxDeltaT();

        // virtual scalar getResidual()
        // {
        //     return residual_;
        // };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} //End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
