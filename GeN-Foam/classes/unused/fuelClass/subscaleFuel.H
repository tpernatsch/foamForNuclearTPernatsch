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
    Foam::subscaleFuel

Description
    Base class used to calculate behavior of fuel. Details of fuel are specified in subclasses.
    It uses a dictionary called fuelProperties

Reference publications
    Carlo Fiorina, Ivor Clifford, Manuele Aufiero, Konstantin Mikityuk,
    "GeN-Foam: a novel OpenFOAM® based multi-physics solver for 2D/3D transient
    analysis of nuclear reactors", Nuclear Engineering and Design, submitted

    Carlo Fiorina, Konstantin Mikityuk, " Application of the new GeN-Foam multi-physics
    solver to the European Sodium Fast Reactor and verification against available codes",
    Proceedings of ICAPP 2015, May 03-06, 2015 - Nice (France), Paper 15226


Author
    Carlo Fiorina <carlo.fiorina@outlook.com; carlo.fiorina@epfl.ch;>

SourceFiles
    subscaleFuel.C

\*---------------------------------------------------------------------------*/

#ifndef subscaleFuel_H
#define subscaleFuel_H

#include "regIOobject.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "fvMatrices.H"
#include "autoPtr.H"
#include "IOFieldFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class subscaleFuel Declaration
\*---------------------------------------------------------------------------*/

class subscaleFuel
:
    public IOdictionary
{
protected:

    // Protected data

        //- Reference to mesh
        const fvMesh& mesh_;

         //- number of fuel zones
        const label fuelZoneNumber_;

        //- list of fuel zones
        mutable PtrList<label> fuelZoneID_;

        //- Solution fields - fuel
        mutable PtrList<IOFieldField<Field, scalar> > Tf_;

        //- Solution fields - clad
        mutable PtrList<IOFieldField<Field, scalar> > Tc_;

        //- submesh size (number of nodes)
        mutable List<label> fuelSubMeshSize_;
        mutable List<label> cladSubMeshSize_;



private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        subscaleFuel
        (
            const subscaleFuel&
        );

        //- Disallow default bitwise assignment
        void operator=(const subscaleFuel&);

public:

    //- Runtime type information
    TypeName("subscaleFuel");

#ifndef SWIG
        declareRunTimeSelectionTable
        (
            autoPtr,
            subscaleFuel,
            dictionary,
            (
                const fvMesh& mesh
            ),
            (mesh)
        );
#endif

    // Constructors

        // Construct from mesh
        subscaleFuel(const fvMesh& mesh);

        //- Selector: Return a reference to the selected model
        static autoPtr<subscaleFuel> New
        (
            const fvMesh& mesh
        );

    // Destructor

        ~subscaleFuel();

    // Member functions

        // Solve for fuel heat sources in fluid given the fluid temperature (Tfl), the heat transfer coefficient between fluid and fuel, and heat source in fuel (qs)
        virtual tmp<volScalarField> heatSources(const volScalarField& Tfl, const volScalarField& hfl, const volScalarField& qs) = 0;

        // Return centerline fuel temperature
        virtual tmp<volScalarField> fuelIn() = 0;

        // Return average fuel temperature
        virtual void fuelAverage(volScalarField& TavFuel) = 0;

        // Return outer radius fuel temperature
        virtual tmp<volScalarField> fuelOut() = 0;

        // Return inner radius clad temperature
        virtual tmp<volScalarField> cladIn() = 0;

        // Return average fuel temperature
        virtual void cladAverage(volScalarField& TavClad) = 0;

        // Return outer radius fuel temperature
        virtual tmp<volScalarField> cladOut() = 0;


};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
