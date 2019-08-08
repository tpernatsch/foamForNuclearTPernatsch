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
    Foam::thermoMechanics

Description
    Base class for neutronic models.

SourceFiles
    thermoMechanics.C

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

#ifndef thermoMechanics_H
#define thermoMechanics_H

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
                          Class thermoMechanics Declaration
\*---------------------------------------------------------------------------*/

class thermoMechanics
:
    public IOdictionary
{

protected:

    // Protected data

        //- Reference to mesh
        fvMesh& mesh_;

		//- Initialize residual
        double initialResidual_;

        //- Mesh displacement
        volVectorField meshDisp_;        

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        thermoMechanics(const thermoMechanics&);

        //- Disallow default bitwise assignment
        void operator=(const thermoMechanics&);


public:

    //- Runtime type information
    TypeName("thermoMechanics")

#ifndef SWIG
        declareRunTimeSelectionTable
        (
            autoPtr,
            thermoMechanics,
            dictionary,
            (
                fvMesh& mesh
            ),
            (mesh)
        );
#endif

    // Constructors

        //- Construct from mesh
        thermoMechanics(fvMesh& mesh);


    // Selectors

        //- Return a reference to the selected model
        static autoPtr<thermoMechanics> New
        (
            fvMesh& mesh
        );


    //- Destructor
    	virtual ~thermoMechanics();


    // Member Functions

        //- Return mesh displacement
        inline const volVectorField& meshDisp() const
        {
            return meshDisp_;
        }    	

        //- Correct/update the properties
        virtual void correct(scalar& residual) = 0;  

        //- Take from another mesh the temperature field
        virtual void getFields(
            const volScalarField& T, 
            const volScalarField& TavFuel,
            const volScalarField& TavClad,
            const meshToMesh& TMToFluid) = 0;  

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
