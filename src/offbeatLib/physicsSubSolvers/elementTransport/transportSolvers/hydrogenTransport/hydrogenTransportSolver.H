/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::transportSolver

Description
    Mother class for the hydrogen transport solvers

SourceFiles
    hydrogenTransportSolver.C

\mainauthor
    M. Reymond (EPFL) 

\contribution
    E.Brunetto, A. Scolaro, C. Fiorina - EPFL
    I. Clifford - PSI

\date 
    October 2022

\*---------------------------------------------------------------------------*/

#ifndef hydrogenTransportSolver_H
#define hydrogenTransportSolver_H

#include "fvMesh.H"
#include "volFields.H"
#include "materials.H"
#include "regIOobject.H"
#include "IOdictionary.H"
#include "globalOptions.H"
#include "transportSolver.H"

#include "TSSpModel.H"
#include "TSSdModel.H"
#include "HDiffusionCoefficientModel.H"
#include "hydrideReorientationModel.H"
#include "physicoChemicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class transportSolver Declaration
\*---------------------------------------------------------------------------*/

class hydrogenTransportSolver
:
    public regIOobject,
    public transportSolver
{
    // Private data
    
    // Private Member Functions
    
        //- Disallow default bitwise copy construct
        hydrogenTransportSolver(const hydrogenTransportSolver&);

        //- Disallow default bitwise assignment
        void operator=(const hydrogenTransportSolver&);
        
protected:
    
    // Protected data

        //- Hydrogen concentration in solid solution (wt.ppm)
        volScalarField Css_;
        
        //- Hydrogen precipitate concentration (wt.ppm)
        volScalarField Cpp_;
        
        //- Total hydrogen concentration (wt.ppm)
        volScalarField Ctot_;
    
        //- Model for terminal solid solubility for dissolution
        // in each cellZone
        PtrList<TSSpModel> TSSpModel_;
        
        //- Model for terminal solid solubility for precipitation
        // in each cellZone
        PtrList<TSSdModel> TSSdModel_;
        
        //- Model to use for diffusion coefficient of hydrogen
        autoPtr<HDiffusionCoefficientModel> HDiffusionCoefficientModel_;

        //- Model to use for hydride reorientation
        autoPtr<hydrideReorientationModel> hydrideReorientationModel_;

        //- Reference to temperature field
        const volScalarField& T_;

        //- Reference to temperature gradient field
        const volVectorField& gradT_;

        //- Hydrogen diffusion coefficient
        volScalarField diffH_;

        //- Hydrogen transport as a function of temperature gradient
        volScalarField soretTerm_;

        //- Threshold concentration for hydrogen precipitation (wt.ppm)
        volScalarField TSSd_;

        //- Threshold concentration for hydrogen dissolution  (wt.ppm)
        volScalarField TSSp_;

        //- Volume fraction of hydrides
        volScalarField hydrideVolFrac_;

        //- Amount of radial hydrides (wt. ppm)
        //- Only available if using cylindrical coordinates
        autoPtr<volScalarField> radialHydrides_;

        //- Fraction of radial hydrides
        autoPtr<volScalarField> fracRHy_;

        //- Molar heat of transport [J/mol]
        dimensionedScalar Q_;

        //- Molar mass of hydrogen
        scalar MH_;
        
        //- Molar mass of zircaloy
        scalar MZr_;

    //- Protected functions
        
        //- Initialise the hydrogen concentration fields
        virtual void initConcentrations();
        
        //- Calculate and report total hydrogen content
        void computeCtot();

        //- Compute atomic fraction of H at the (delta / (alpha + delta)) boundary
        virtual scalar hydrideBoundaryFraction(const scalar&);

        //- Compute hydride volume fraction
        virtual scalar hydrideVolFrac(const scalar&, const scalar&, const scalar&);
  
public: 

    //- Runtime type information
    TypeName("hydrogenTransport");  
        
    // Constructors
        //- Construct from mesh, materials and dictionary
        hydrogenTransportSolver
        (
            fvMesh& mesh,
            const materials& mat,
            const dictionary& elementTransportDict,
            const word& solverName
        ); 

    //- Destructor
    virtual ~hydrogenTransportSolver();

    // Access functions
    
        //- Name of primary solution field
        virtual const word& fieldName()
        {
            return Css_.name();
        }
        
        //- Return the molar heat of transport [J/mol]
        inline const dimensionedScalar& Q() const
        {
            return Q_;
        }
        
        //- Return the temperature field
        inline const volScalarField& T() const
        {
            return T_;
        }
        
        //- Return the diffusion coefficient field
        inline const volScalarField& diffH() const
        {
            return diffH_;
        }

    // Member Functions
     
        //- Calculate the updated hydrogen concentrations
        virtual void correct();

        //- Next time step size according to solver criterion
        virtual scalar nextDeltaT();
        
        //- Virtual function needed by regIOobject
        virtual bool writeData (Ostream& os) const
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
