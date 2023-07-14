/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2306                                                  |
|    Copyright 2011-2016 OpenFOAM Foundation, 2017-2022 OpenCFD Ltd.         |
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
    - \f$ \phi_g \f$ is the neutron flux for the gth energy group
    - \f$ C_k \f$ is the concentration of the kth delayed neutron precursor
        group
    - \f$ D_g \f$ is the neutron diffusion coefficient for the gth energy group
    - \f$ \chi_{p,g} \f$ is the prompt neutron yield for the gth energy group
    - \f$ \chi_{d,g} \f$ is the delayed neutron yield for the ith energy group
    - \f$ k_{eff} \f$ is the effective multiplication factor
    - \f$ \beta_{tot} \f$ is the total delayed neutron fraction
    - \f$ \beta_k \f$ is the delayed neutron fraction for the kth delayed 
        neutron precursor group
    - \f$ \nu \f$ is the average number of neutrons per fission
    - \f$ \Sigma_{f,g} \f$ is the fission cross section for the gth energy group
    - \f$ \Sigma_{r,g} \f$ is the removal (disappearance) cross section for the 
        gth energy group
    - \f$ \Sigma_{h \rightarrow g} \f$ is the group-transfer cross section from 
        the hth to the gth energy group
    - \f$ S_{n,g} \f$ is the fission neutron source from neutron energy groups
        others than the gth
    - \f$ S_d \f$ is the delayed neutron source
    - \f$ S_{s,g}  \f$ is the scattering neutron source from neutron energy 
        groups others than the gth
    - \f$ S_g \f$ is the external neutron source for the gth energy group
    - \f$ \boldsymbol{u} \f$ is the fuel velocity
    - \f$v_g \f$ is the average neutron velocity for the ith energy group
    - \f$ \lambda_k \f$ is the decay constant for the kth delayed neutron 
        precursor group


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

    Reference publications

    Carlo Fiorina, Ivor Clifford, Manuele Aufiero, Konstantin Mikityuk,
    "GeN-Foam: a novel OpenFOAM® based multi-physics solver for 2D/3D transient
    analysis of nuclear reactors", Nuclear Engineering and Design, submitted. 
    
    Carlo Fiorina, Konstantin Mikityuk, " Application of the new GeN-Foam 
    multi-physics solver to the European Sodium Fast Reactor and verification 
    against available codes", Proceedings of ICAPP 2015, May 03-06, 2015 -
    Nice (France), Paper 15226

    Carlo Fiorina, Nordine Kerkar, Konstantin Mikityuk, Pablo Rubiolo, Andreas 
    Pautz, "Development and verification of the neutron diffusion solver for the 
    GeN-Foam multi-physics platform", Annals of Nuclear Energy, 96, 2016, 
    pages 212-222, DOI: http://dx.doi.org/10.1016/j.anucene.2016.05.023.

    Author of this file (and associated .C or included .H files): 
    Carlo Fiorina <carlo.fiorina@outlook.com; carlo.fiorina@epfl.ch;>, EPFL (Switzerland)   
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
    volScalarField oneGroupFlux_;
    volScalarField neutroSource_;
    volScalarField delayedNeutroSource_;
    volScalarField scatteringSourceExtra_;

    //- Ptrs to coupling fields (i.e. defined on their original mesh)
    const volScalarField* TFuelOrig_;
    const volScalarField* TCladOrig_;
    const volScalarField* TCoolOrig_;
    const volScalarField* rhoCoolOrig_;
    const volVectorField* UOrig_;
    const volScalarField* alphaOrig_; //- Total volume fraction of fluid(s)
    const volScalarField* alphatOrig_; //- Thermal diffusivity in kg/m/s
    const volScalarField* muOrig_; //- Molecular viscosity in kg/m/s (== Pa*s)

    //- Mapped fields
    volScalarField TFuel_;
    volScalarField TClad_;
    volScalarField TCool_;
    volScalarField rhoCool_;
    
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
        virtual const volScalarField& oneGroupFlux() const
        {
            return oneGroupFlux_;
        }

        //- Return power
        virtual scalar power() const
        {
            scalar pTot = 0;

            forAll(flux_,energyI)
            {
                pTot += fvc::domainIntegrate
                (
                    flux_[energyI]*xs_.sigmaPow()[energyI]
                ).value();
            }
            return pTot;
        }

        virtual void getCouplingFieldRefs
        (
            const objectRegistry& src,
            const meshToMesh& neutroToFluid
        );   

        virtual void interpolateCouplingFields
        (
            const meshToMesh& neutroToFluid
        );

        //- Correct/update the properties
        virtual void correct
        (
            scalar& residual,
            label couplingIter  
        );  
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
