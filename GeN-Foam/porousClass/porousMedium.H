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
    Foam::porousMedium

Description
    Base class for porous media models.

SourceFiles
    porousMedium.C

Reference publications
	Carlo Fiorina, Ivor Clifford, Manuele Aufiero, Konstantin Mikityuk,
	"GeN-Foam: a novel OpenFOAM® based multi-physics solver for 2D/3D transient
	analysis of nuclear reactors", Nuclear Engineering and Design, submitted

	Carlo Fiorina, Konstantin Mikityuk, " Application of the new GeN-Foam multi-physics
	solver to the European Sodium Fast Reactor and verification against available codes",
	Proceedings of ICAPP 2015, May 03-06, 2015 - Nice (France), Paper 15226


Author
    Carlo Fiorina <carlo.fiorina@hotmail.it; carlo.fiorina@psi.ch; carlo.fiorina@polimi.it; 
                   carlo.fiorina@epfl.ch;>

\*---------------------------------------------------------------------------*/

#ifndef porousMedium_H
#define porousMedium_H

#include "volFields.H"
#include "rhoThermo.H"
#include "fvMatrices.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "coordinateSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class porousMedium Declaration
\*---------------------------------------------------------------------------*/

class porousMedium
:
    public IOdictionary //eredita tutto cio' che c'e' in IOdictionary
{

protected:

    // Protected data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Reference to thermophysical properties
		const rhoThermo& thermo_;

        //- Reference to density field
        const volScalarField& rho_;

        //- Reference to velocity field
        const volVectorField& U_;

        //- Reference to viscosity field
        const volScalarField& mu_;

        //- Void Fraction
        volScalarField gamma_;

        //- Volumetric Area of subscale structures
        volScalarField volumetricAreaSS_;

        //- Volumetric Area of fuel
        volScalarField volumetricAreaFuel_;

        //- rotation tensor from local to global coordinate system
        volTensorField rotate_;

        //- heat transfer coeff from the porous structure to some external environment
        volScalarField externalHeatTransferCoefficient_;

        //- temperature of the external environment
        volScalarField externalT_;

        //- heat source coming from somewhere outside the porous medium
        volScalarField externalVolHeatSource_;

        //- rho*cp of the pororus medium (useful for transients)
        volScalarField externalRhoCp_;

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        porousMedium(const porousMedium&);

        //- Disallow default bitwise assignment
        void operator=(const porousMedium&);


public:

    //- Runtime type information
    TypeName("porousMedium")

#ifndef SWIG
        declareRunTimeSelectionTable
        (
            autoPtr,
            porousMedium,
            dictionary,
            (
                const volScalarField& rho,
                const volVectorField& U,
				const rhoThermo& thermo
            ),
            (rho, U, thermo)
        );
#endif

    // Constructors

        //- Construct from mesh
        porousMedium
        (
            const volScalarField& rho, // in principle not necessary
            const volVectorField& U,
	    const rhoThermo& thermo
        );


    // Selectors

        //- Return a reference to the selected model
        static autoPtr<porousMedium> New
        (
            const volScalarField& rho,
            const volVectorField& U,
	    	const rhoThermo& thermo
        );


    //- Destructor
    virtual ~porousMedium();


    // Member Functions

        // Access

            //- Return the void fraction
            inline const volScalarField& voidFraction() const //the last const says that this function will not change anything in the class
            {
                return this->gamma_;
            }

            //- Return the Volumetric Area of subscale structures
            inline const volScalarField& volumetricAreaSS() const //the last const says that this function will not change anything in the class
            {
                return this->volumetricAreaSS_;
            }

            //- Return the Volumetric Area of fuel
            inline const volScalarField& volumetricAreaFuel() const 
            {
                return this->volumetricAreaFuel_;
            }

            //- Return rotate (global to local coord system)
            inline const volTensorField& rotate() const
            {
                return rotate_;
            }

            //- Return externalHeatTransferCoefficient
            inline const volScalarField& externalHeatTransferCoefficient() const
            {
                return externalHeatTransferCoefficient_;
            }

            //- Return externalT
            inline const volScalarField& externalT() const
            {
                return externalT_;
            }

            //- Return externalVolHeatSource
            inline const volScalarField& externalVolHeatSource() const
            {
                return externalVolHeatSource_;
            }

            //- Return externalVolHeatSource
            inline const volScalarField& externalRhoCp() const
            {
                return externalRhoCp_;
            }

            //- Return the drag coefficient
            virtual tmp<volTensorField> dragCoeff() const = 0;   // use virtual to make him use the most derived version of this function. =0 makes it a pure virtual function, i.e. you must fefine it in a derived class while you are allowed not to define it here

            //- Return the pumpMomentumSource
            virtual tmp<volVectorField> pump() const = 0;

            //- Return the Reynolds number
            virtual tmp<volScalarField> ReynoldsNumber() const = 0;

            //- Return the heat transfer coefficient
            virtual tmp<volScalarField> heatTransferCoefficient() const = 0;

            //- Return the explicit drag force
            tmp<volVectorField> explicitMomentumSource() const;

            //- Return the implicit drag force
            fvVectorMatrix semiImplicitMomentumSource(volVectorField& U) const;


        // Edit

            //- Update coefficients
            virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
