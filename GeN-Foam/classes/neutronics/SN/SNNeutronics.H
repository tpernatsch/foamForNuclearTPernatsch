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
    Foam::neutronics::SNNeutronics

Description
    subclass for diffusion neutronics.

SourceFiles
    SNNeutronics.C

Reference publications
    Carlo Fiorina, Ivor Clifford, Manuele Aufiero, Konstantin Mikityuk,
    "GeN-Foam: a novel OpenFOAM® based multi-physics solver for 2D/3D transient
    analysis of nuclear reactors", Nuclear Engineering and Design, submitted

    Carlo Fiorina, Konstantin Mikityuk, " Application of the new GeN-Foam multi-physics
    solver to the European Sodium Fast Reactor and verification against available codes",
    Proceedings of ICAPP 2015, May 03-06, 2015 - Nice (France), Paper 15226


Author
    Derived by Carlo Fiorina <carlo.fiorina@outlook.com; carlo.fiorina@epfl.ch;>
    from the work of Manuele Aufiero, "D EVELOPMENT OF ADVANCED SIMULATION TOOLS FOR
    CIRCULATING - FUEL NUCLEAR REACTORS", PhD Dissertation, Politecnico Di Milano


\*---------------------------------------------------------------------------*/

#ifndef SNNeutronics_H
#define SNNeutronics_H

#include "neutronics.H"
#include "XSLowMem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class SNNeutronics Declaration
\*---------------------------------------------------------------------------*/

class SNNeutronics
:
    public neutronics
{
protected:

    // Cross sections
    XSLowMem xs_;

    //quadrature Set
    IOdictionary quadratureSet_;
    label discreteDirections_;
    vectorField directionVersors_;
    scalarField directionWeights_;
    PtrList<scalarSquareMatrix> legendreMatrices_;
    surfaceScalarField facePhi_;
    PtrList<surfaceScalarField>  facePhis_;

    //Fields
    List<PtrList<volScalarField> > angularFlux_;
    PtrList<volScalarField> prec_;//prec per total volume
    PtrList<volScalarField> precStar_;//prec per volume salt
    volScalarField defaultFlux_;
    volScalarField defaultPrec_;

    //Derived fields
    PtrList<volScalarField> flux_;
    volScalarField oneGroupFlux_;
    volScalarField neutroSource_;
    volScalarField delayedNeutroSource_;

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

        //- Access

            label discreteDirections() const
            {
                return discreteDirections_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
