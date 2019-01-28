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

    Carlo Fiorina, Konstantin Mikityuk, " Application of the new GeN-Foam multi-physics
    solver to the European Sodium Fast Reactor and verification against available codes",
    Proceedings of ICAPP 2015, May 03-06, 2015 - Nice (France), Paper 15226


Author
    Carlo Fiorina <carlo.fiorina@outlook.com; carlo.fiorina@epfl.ch;>

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
#include "fvMatrixExt.H"
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

    // Protected data

        //- Reference to mesh
        fvMesh& mesh_;

        IOdictionary reactorState_;

        //- k effective
		scalar keff_;

        //- Target power
		scalar pTarget_;

		//- Eigenvalue flag
		bool eigenvalueNeutronics_;

		//- Initialize residual
        double initialResidual_;

        //- Vol power
        volScalarField volFuelPower_;

        //- Mesh Displacement
        volVectorField Disp_;

private:

    // Private Member Functions

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

    // Constructors

        //- Construct from mesh
        neutronics(fvMesh& mesh);


    // Selectors

        //- Return a reference to the selected model
        static autoPtr<neutronics> New
        (
            fvMesh& mesh
        );


    //- Destructor
    	virtual ~neutronics();


    // Member Functions

        //- Return externalVolHeatSource
        inline const bool& eigenvalue() const
        {
            return eigenvalueNeutronics_;
        }    	

        //- Return one group flux
        inline const scalar& keff() const
        {
            return keff_;
        }

        //- Return externalVolHeatSource
        inline const volScalarField& volFuelPower() const
        {
            return volFuelPower_;
        }    	

        //- Return mesh displacement
        inline const volVectorField& Disp() const
        {
            return Disp_;
        }    	

        //- Return fluxes
        virtual inline const PtrList<volScalarField> fluxes() const = 0;

        //- Return precursors
        virtual inline const PtrList<volScalarField> precursors() const = 0;

        //- Return one group flux
        virtual inline const volScalarField& Flux() const = 0;

        //- Return power
        virtual scalar power() const = 0;  

        //- Correct/update the properties
        virtual void correct(const label couplingIter, scalar& residual, const bool& liquidFuel ) = 0;  

        //- Take from another mesh the fields for XS parametrization
        virtual void getFields(
            const volScalarField& TfuelOrig, 
            const volScalarField& TcladOrig, 
            const volScalarField& rhoCoolOrig, 
            const volScalarField& TCoolOrig, 
            const meshToMesh& neutroToFluid) = 0;  

        virtual void getFieldsLiquidFuel(
            const volVectorField& UOrig, 
            const volScalarField& porosityOrig, 
            const volScalarField& TfuelOrig, 
            const volScalarField& TcladOrig, 
            const volScalarField& rhoCoolOrig, 
            const volScalarField& TCoolOrig,
            const volScalarField& muOrig, 
            const volScalarField& alphatOrig, 
            const meshToMesh& neutroToFluid) = 0; 

        //- Deform mesh
        virtual void deformMesh(const meshToMesh& TMToNeutro,const volVectorField& DispOrig) = 0; 

        //- Calculate initialResidual_
        inline double initialResidual() const
        {
            return initialResidual_;
        };

        //- Write data to time directory
        virtual bool writeData(Ostream& os) const ;    

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
