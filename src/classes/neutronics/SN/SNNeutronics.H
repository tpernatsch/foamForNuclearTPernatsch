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
    Foam::SNNeutronics

Description
    Subclass for SN neutronics (Discrete Ordinate) [@FIORINA201524] [@Fiorina2015ApplicationCodes].

$$
\begin{align}
    \frac{1}{v_i} \frac{d \Omega_j \cdot \phi_{i,j}}{dt}
        + \nabla \cdot (\Omega_j \cdot \phi_{i,j})
        + \Sigma_{t,i} \phi_{i,j}
    & = S_{i,j}^I + S_{i,j}^A + S_i^D
        \quad \quad \forall i \wedge j \in [1; G] \times [1; M] \\
    \frac{\partial C_k}{\partial t} + \nabla \cdot (\mathbf{u} C_k) & =
        \frac{\beta_k}{k_{eff}} \sum_{i'=1}^G \nu \Sigma_{f,i'} \phi_{i'}
        - \lambda_k C_k
        \quad \forall k \in [1; N]
\end{align}
$$

with the external source:

$$
\begin{align}
    S_{i,j}^I & = \sum_{i'=1}^G \sum_{j'=1}^M \omega_{j'} \phi_{i',j'} \left(
        \frac{1}{k_{eff}} \nu \Sigma_{f,i'} \chi_{p,i} + \Sigma_{s,i' \rightarrow i}
    \right) \\
    S_{i,j}^A & = \sum_{i'=1}^G \sum_{j'=1}^M \omega_{j'} \phi_{i',j'} \Sigma_{s,i' \rightarrow i} \left(
        \sum_{l=1}^p s_{i' \rightarrow i}^l \frac{2l+1}{2} P_l (\Omega_j \cdot \Omega_{j'})
    \right) \\
    S_i^D & = \chi_{d,i} \sum_k \lambda_k C_k
\end{align}
$$

The anisotropic scattering terms are considered adopting a pth order Legendre
polynomial expansion of the scattering cosine \f$ (\Omega_j \cdot \Omega_{j'}) \f$

\vartable
    \phi_{i,j}              | Neutron flux if the energy group i and direction j
    C_k                     | Precursor k
    S_{i,j}^I               | Isotropic source term
    S_{i,j}^A               | Anisotropic source term
    S_i^D                   | Delayed neutron source term
    \Sigma_{t,i}            | Total cross-section of the energy group i
    \nu \Sigma_{f,i'}       | Fission neutron production cross-sections of the energy group i'
    \Sigma_{s,i' \rightarrow i} | Scattering cross-section from energy i' to energy i
    k_{eff}                 | Multiplication factor
    \chi_{p,i}              | Prompt emission neutrons emitted in the energy group i
    \chi_{d,i}              | Delayed emission neutrons emitted in the energy group i
    \lambda_k               | Decay constant of the precursor group k
    \beta_k        | Delayed neutron fraction for the kth delayed neutron precursor group
    \omega_{j'}             | Normalized quadrature weight associated to the direction \f$ \Omega_{j'} \f$
    s_{i \rightarrow i'}^l  | lth Legendre coefficient of the scattering cross-section \f$ \Sigma_{s,i \rightarrow i'} \f$
    P_l                     | lth Legendre polynomial
    \mathbf{u}          | Fuel velocity
    v_i                     | Average neutron velocity of the energy group i
    \Omega_j                | Unit vector of the angular direction j
\endvartable

<b>Reference publications</b>

<ol type="1">
    <li> Carlo Fiorina, Ivor Clifford, Manuele Aufiero, Konstantin Mikityuk,
    GeN-Foam: a novel OpenFOAM® based multi-physics solver for 2D/3D transient
    analysis of nuclear reactors, Nuclear Engineering and Design, Volume 294,
    2015, Pages 24-37, ISSN 0029-5493, https://doi.org/10.1016/j.nucengdes.2015.05.035.

    <li> Carlo Fiorina, Konstantin Mikityuk, "Application of the new GeN-Foam multi-physics
    solver to the European Sodium Fast Reactor and verification against available codes",
    Proceedings of ICAPP 2015, May 03-06, 2015 - Nice (France), Paper 15226
</ol>

\mainauthor of this file (and associated .C or included .H files):
    Derived by Carlo Fiorina <carlo.fiorina@outlook.com; carlo.fiorina@epfl.ch;> (EPFL, Switzerland)
    from the work of Manuele Aufiero, "DEVELOPMENT OF ADVANCED SIMULATION TOOLS FOR
    CIRCULATING - FUEL NUCLEAR REACTORS", PhD Dissertation, Politecnico Di Milano

SourceFiles
    SNNeutronics.C

\*---------------------------------------------------------------------------*/

#ifndef SNNeutronics_H
#define SNNeutronics_H

#include "neutronics.H"
#include "XS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
/*---------------------------------------------------------------------------*\
                          Class SNNeutronics Declaration
\*---------------------------------------------------------------------------*/

class SNNeutronics
:
    public neutronics
{
protected:

    //- Cross sections
    XS xs_;

    // Quadrature Set
    IOdictionary quadratureSet_;
    label discreteDirections_;
    vectorField directionVersors_;
    scalarField directionWeights_;
    PtrList<scalarSquareMatrix> legendreMatrices_;
    surfaceScalarField facePhi_;
    PtrList<surfaceScalarField> facePhis_;

    //- Fields
    List<PtrList<volScalarField>> angularFlux_;

    //- Precursors per total volume
    PtrList<volScalarField> prec_;

    //- Precursors per volume salt
    PtrList<volScalarField> precStar_;
    volScalarField defaultFlux_;
    volScalarField defaultPrec_;

    //Derived fields
    PtrList<volScalarField> flux_;
    volScalarField neutroSource_;
    volScalarField delayedNeutroSource_;

    //- For Aitken eigenvalue
    scalar keff0_;
    scalar keff1_;
    scalar keff2_;
    int aitkenIterCounter_;
    int aitkenIterNo_;

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

    // Private data

    // Private Member Functions

        //- Disallow default bitwise copy construct
        SNNeutronics(const SNNeutronics&);

        //- Disallow default bitwise assignment
        void operator=(const SNNeutronics&);


public:

    //- Runtime type information
        TypeName("SNNeutronics");


    // Constructors

        //- Construct from mesh, read data from IOdictionaries
        SNNeutronics(fvMesh& mesh);

    //- Destructor
        virtual ~SNNeutronics();


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

        //- Access

            label discreteDirections() const
            {
                return discreteDirections_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} //End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
