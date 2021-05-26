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
    Foam::neutronics

Description
    Base class for neutronic models.

SourceFiles
    neutronics.C

Reference publications
    Carlo Fiorina, Ivor Clifford, Manuele Aufiero, Konstantin Mikityuk,
    "GeN-Foam: a novel OpenFOAM® based multi-physics solver for 2D/3D transient
    analysis of nuclear reactors", Nuclear Engineering and Design, submitted

    Carlo Fiorina, Konstantin Mikityuk, " Application of the new GeN-Foam 
    multi-physics solver to the European Sodium Fast Reactor and verification 
    against available codes", Proceedings of ICAPP 2015, May 03-06, 2015 -
    Nice (France), Paper 15226

Author of this file (and associated .C or included .H files)
    Carlo Fiorina <carlo.fiorina@outlook.com; carlo.fiorina@epfl.ch;>, EPFL (Swizerland)
    

\*---------------------------------------------------------------------------*/

#ifndef neutronics_H
#define neutronics_H

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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class neutronics Declaration
\*---------------------------------------------------------------------------*/

class neutronics
:
    public IOdictionary
{

protected:

    //- Reference to mesh
    fvMesh& mesh_;

    //- Dictionary from which keff, pTarget are read
    IOdictionary reactorState_;

    //- k effective
    scalar keff_;

    //- Target power
    scalar pTarget_;

    //- Fuel power density
    volScalarField powerDensity_;

    //- Mesh displacement
    volVectorField disp_;

    //- Initialize residual
    scalar initialResidual_;

    //- Eigenvalue flag
    bool eigenvalueNeutronics_;

    //- Liquid fuel (e.g. true for MSRs)
    const bool liquidFuel_;

private:

    //- Disallow default bitwise copy construct
    neutronics(const neutronics&);

    //- Disallow default bitwise assignment
    void operator=(const neutronics&);

public:

    //- Runtime type information
    TypeName("neutronics")

#ifndef SWIG
        declareRunTimeSelectionTable
        (
            autoPtr,
            neutronics,
            dictionary,
            (
                fvMesh& mesh
            ),
            (mesh)
        );
#endif

    //- Constructors

        //- Construct from mesh
        neutronics(fvMesh& mesh);

    //- Selectors

        //- Return a reference to the selected model
        static autoPtr<neutronics> New
        (
            fvMesh& mesh
        );

    //- Destructor
        virtual ~neutronics(){}

    //- Member Functions

        //- Deform mesh
        virtual void deformMesh
        (
            const meshToMesh& TMToNeutro,
            const volVectorField& dispOrig
        ); 

        //- Return fluxes
        virtual const PtrList<volScalarField> fluxes() const = 0;

        //- Return precursors
        virtual const PtrList<volScalarField> precursors() const = 0;

        //- Return one group flux
        virtual const volScalarField& oneGroupFlux() const = 0;

        //- Return power
        virtual scalar power() const = 0;  

        //- Substitutes the old getFields. Set local class references to global
        //  coupling fields, to be called only once after construction. The 
        //  fields are set via lookup in the src objectRegistry (i.e. 
        //  fluidMeshNB for the current state of GeN-Foam)
        virtual void getCouplingFieldRefs
        (
            const objectRegistry& src,
            const meshToMesh& neutroToFluid
        ) = 0;

        //- The ond getFields, yet operates on the local refs without having
        //  the need for all the fields to be passed. This effectively bypasses
        //  the old limitation of getFields, i.e. if someone wanted to
        //  implement a new neutronics model that requires more fields for the
        //  coupling (as it happend with pointKinetics), then they'd need to
        //  modify the getFields function in ALL of the other neutronics model
        //  too. Splitting getFields into a getCouplingFieldRefs and 
        //  interpolateCouplingFields allows one to be as flexible as possible
        //  with new classes implementations
        virtual void interpolateCouplingFields
        (
            const meshToMesh& neutroToFluid
        ) = 0;

        //- Correct/update the properties
        virtual void correct
        (
            scalar& residual,
            label couplingIter
        ) = 0;  

        //- Required to remove "abstractedness" from IOdictionary
        virtual bool writeData(Ostream& os) const
        {
            return os.good();
        }

        //- Access

            scalar keff() const
            {
                return keff_;
            }

            const volScalarField& powerDensity() const
            {
                return powerDensity_;
            }       

            const volVectorField& disp() const
            {
                return disp_;
            }       

            bool eigenvalue() const
            {
                return eigenvalueNeutronics_;
            } 

            scalar initialResidual() const
            {
                return initialResidual_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
