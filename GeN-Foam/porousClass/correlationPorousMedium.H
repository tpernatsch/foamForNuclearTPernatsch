/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::porousMedium::correlationPorousMedium

Description
    subclass of porousMedium. Created to evaluate pressure drops and
    heat transfer in a porous medium using available correlations.
    Different correlations can be used for three principal directions
    of a local coordinate system defined by the user. The possibility
    also exists to use different correlations for laminar and turbulent
    flow.

SourceFiles
    correlationPorousMedium.C

Reference publications
    Carlo Fiorina, Ivor Clifford, Manuele Aufiero, Konstantin Mikityuk,
    "GeN-Foam: a novel OpenFOAM® based multi-physics solver for 2D/3D transient
    analysis of nuclear reactors", Nuclear Engineering and Design, submitted

    Carlo Fiorina, Konstantin Mikityuk, " Application of the new GeN-Foam multi-physics
    solver to the European Sodium Fast Reactor and verification against available codes",
    Proceedings of ICAPP 2015, May 03-06, 2015 - Nice (France), Paper 15226


Author
    Carlo Fiorina <carlo.fiorina@outlook.com; carlo.fiorina@epfl.ch;>

\*---------------------------------------------------------------------------*/

#ifndef correlationPorousMedium_H
#define correlationPorousMedium_H

#include "porousMedium.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class correlationPorousMedium Declaration
\*---------------------------------------------------------------------------*/

class correlationPorousMedium
:
    public porousMedium//class inheritance
{
protected:


        //- Transition Re numbers
        volVectorField reynoldsTurb_;
        volVectorField reynoldsLam_;

        //- A in (A*Re^a) for Darcy friction
        volVectorField darcyConstTurb_;
        volVectorField darcyConstLam_;

        //- a in (A*Re^a) for Darcy friction
        volVectorField darcyExpTurb_;
        volVectorField darcyExpLam_;

        //- A in (A*Re^a*Pr^b +B) for nusselt
        volVectorField nusseltConstTurb1_;
        volVectorField nusseltConstLam1_;

        //- B in (A*Re^a*Pr^b +B) for nusselt
        volVectorField nusseltConstTurb2_;
        volVectorField nusseltConstLam2_;

        //- a in (A*Re^a*Pr^b +B) for nusselt
        volVectorField nusseltReExpTurb_;
        volVectorField nusseltReExpLam_;

        //- b in (A*Re^a*Pr^b +B) for nusselt
        volVectorField nusseltPrExpTurb_;
        volVectorField nusseltPrExpLam_;

        //- hydraulicDiameter
        volScalarField hydraulicDiameter_;

        //- hydraulicDiameter
        volScalarField hydraulicDiameterStructure_;

        //- pumpMomentumSource
        volVectorField pumpMomentumSource_;

        //- A in (A*Re^a) for turbulent intensity
        volVectorField turbulenceIntensityConst_;

        //- a in (A*Re^a) for turbulent intensity
        volVectorField turbulenceIntensityExp_;

        //- A in (A*Dh) for turbulent length scale
        volVectorField turbulenceLengthScaleConst_;

        //- k and epsilon convergence rate (1/cm) to their asymptotic value in the porous medium
        volVectorField kepsilonConvergenceRate_;


private:

    // Private data

    // Private Member Functions

        //- Disallow default bitwise copy construct
        correlationPorousMedium(const correlationPorousMedium&);

        //- Disallow default bitwise assignment
        void operator=(const correlationPorousMedium&);


public:

    //- Runtime type information
    TypeName("correlation");


    // Constructors

        //- Construct from mesh, read coefficients from IOdictionary
        correlationPorousMedium
        (
            const volScalarField& rho,
            const volVectorField& U,
            const rhoThermo& thermo
        );

        //- Construct from coefficients
        correlationPorousMedium
        (
            const volScalarField& rho,
            const volVectorField& U,
            const rhoThermo& thermo,
            const dimensionedScalar& gamma,
            const dimensionedScalar& volumetricAreaSS,
            const dimensionedScalar& volumetricAreaFuel,
            const dimensionedScalar& volumetricAreaHX,
            const dimensionedTensor& rotate,
            const dimensionedScalar& externalHeatTransferCoefficient,
            const dimensionedScalar& externalT,
            const dimensionedScalar& externalVolHeatSource,
            const dimensionedScalar& externalRhoCp,
            const dimensionedVector& reynoldsTurb,
            const dimensionedVector& reynoldsLam,
            const dimensionedVector& darcyConstTurb,
            const dimensionedVector& darcyConstLam,
            const dimensionedVector& darcyExpTurb,
            const dimensionedVector& darcyExpLam,
            const dimensionedVector& nusseltConstTurb1,
            const dimensionedVector& nusseltConstLam1,
            const dimensionedVector& nusseltConstTurb2,
            const dimensionedVector& nusseltConstLam2,
            const dimensionedVector& nusseltReExpTurb,
            const dimensionedVector& nusseltReExpLam,
            const dimensionedVector& nusseltPrExpTurb,
            const dimensionedVector& nusseltPrExpLam,
            const dimensionedScalar& hydraulicDiameter,
            const dimensionedScalar& hydraulicDiameterStructure,
            const dimensionedVector& pumpMomentumSource,
            const dimensionedVector& turbulenceIntensityConst_,
            const dimensionedVector& turbulenceIntensityExp_,
            const dimensionedVector& turbulenceLengthScaleConst_,
            const dimensionedVector& kepsilonConvergenceRate_
        );

    //- Destructor
    virtual ~correlationPorousMedium();


    // Member Functions

        // Access

            //- Return turbulent transition Re
            inline const volVectorField& reynoldsTurb() const
            {
                return reynoldsTurb_;
            }

            //- Return lamina transition Re
            inline const volVectorField& reynoldsLam() const
            {
                return reynoldsLam_;
            }

            //- Return darcyConst
            inline const volVectorField& darcyConstTurb() const
            {
                return darcyConstTurb_;
            }

            //- Return darcyConst
            inline const volVectorField& darcyConstLam() const
            {
                return darcyConstLam_;
            }

            //- Return darcyExp
            inline const volVectorField& darcyExpTurb() const
            {
                return darcyExpTurb_;
            }

            //- Return darcyExp
            inline const volVectorField& darcyExpLam() const
            {
                return darcyExpLam_;
            }

            //- Return nusseltConst1
            inline const volVectorField& nusseltConstTurb1() const
            {
                return nusseltConstTurb1_;
            }

            //- Return nusseltConst1
            inline const volVectorField& nusseltConstLam1() const
            {
                return nusseltConstLam1_;
            }

            //- Return nusseltConst2
            inline const volVectorField& nusseltConstTurb2() const
            {
                return nusseltConstTurb2_;
            }

            //- Return nusseltConst2
            inline const volVectorField& nusseltConstLam2() const
            {
                return nusseltConstLam2_;
            }

            //- Return nusseltReExpTurb
            inline const volVectorField& nusseltReExpTurb() const
            {
                return nusseltReExpTurb_;
            }

            //- Return nusseltReExpLam
            inline const volVectorField& nusseltReExpLam() const
            {
                return nusseltReExpLam_;
            }

            //- Return nusseltPrExpTurb
            inline const volVectorField& nusseltPrExpTurb() const
            {
                return nusseltPrExpTurb_;
            }

            //- Return nusseltPrExpLam
            inline const volVectorField& nusseltPrExpLam() const
            {
                return nusseltPrExpLam_;
            }

            //- Return hydraulicDiameter
            inline const volScalarField& hydraulicDiameter() const
            {
                return hydraulicDiameter_;
            }

            //- Return hydraulicDiameter
            inline const volScalarField& hydraulicDiameterStructure() const
            {
                return hydraulicDiameterStructure_;
            }

            //- Return pumpMomentumSource
            inline const volVectorField& pumpMomentumSource() const
            {
                return pumpMomentumSource_;
            }

            //- Return turbulenceIntensityConst
            inline const volVectorField& turbulenceIntensityConst() const
            {
                return turbulenceIntensityConst_;
            }

            //- Return turbulenceIntensityExp
            inline const volVectorField& turbulenceIntensityExp() const
            {
                return turbulenceIntensityExp_;
            }

            //- Return turbulenceLengthScaleConst
            inline const volVectorField& turbulenceLengthScaleConst() const
            {
                return turbulenceLengthScaleConst_;
            }

        // Calc

            //-  darcy friction coefficient
            tmp<volVectorField> darcyFrictionFactor() const;

            //-  equilibrium k
            tmp<volScalarField> equilibriumK() const;

            //-  equilibrium epsilon
            tmp<volScalarField> equilibriumEpsilon() const;

            //- Return kepsilonConvergenceRate, in 1/s
            tmp<volScalarField> kepsilonConvergenceRate() const ;


        // From porous medium

            //- Return the drag coefficient tensor
            virtual tmp<volTensorField> dragCoeff() const ;

            //- Return the pumpMomentumSource
            virtual tmp<volVectorField> pump() const ;

            //- Prandtl number
            virtual tmp<volScalarField>  PrandtlNumber() const ;

            //- Reynolds number
            virtual tmp<volScalarField>  ReynoldsNumber() const ;

            //- Return the heat transfer coefficient
            virtual tmp<volScalarField> heatTransferCoefficient() const ;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
