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
    Foam::gapGasModel

Description
	As the name implies, the `none` gapGas class in OFFBEAT neglects the presence
    of the gap gas.

    \warning        
    Some OFFBEAT custom boundary conditions for gap heat trasnfer or contact assume 
    the presence of a gap gas model, and will prompt a Fatal Error message in case 
    the `none` gapGas model is selected.

    The class derives from regIOobject. It writes/reads a file named `gapGas`
    at each time step.

________________________________________________________________________________

Usage
	Here is a code snippet of the `solverDict` to be used for activating the
	`none` gapGas model:

    \verbatim
        gapGas none;
    \endverbatim

SourceFiles
    gapGasModel.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef gapGasModel_H
#define gapGasModel_H

#include "Time.H"
#include "IOdictionary.H"
#include "fvMesh.H"

#include "volFields.H"
#include "surfaceFields.H"
#include "Switch.H"

#include "materials.H"
#include "fissionGasRelease.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class gapGasModel Declaration
\*---------------------------------------------------------------------------*/

class gapGasModel
:
    public regIOobject
{
// Private data
            
    //- Disallow default bitwise copy construct
    gapGasModel(const gapGasModel&);

    //- Disallow default bitwise assignment
    void operator=(const gapGasModel&);
        
protected:
    
    //- Reference to mesh
    const fvMesh& mesh_;

    //- Reference to the materials class
    const materials& mat_;
    
    //- Reference to temperature field
    const volScalarField& T_;

    //- Reference to displacement (or incremental displacement) field
    const volVectorField& DorDD_;

    //- Reference to the fgr class
    const fissionGasRelease& fgr_;

    //- gapGas dictionary
    const dictionary gapGasDict_;
    
    //- Mixture component names
    wordList species_;  
    
    //- Species molecolar weights
    scalarList speciesW_;
    
    //- Gap volume [m3]
    scalar gapV_; 
    
    //- Gap temperature [K]
    scalar gapT_;         
    
    //- Gas pressure [Pa]
    scalar gasP_;     

    //- Gas mass [kg]
    scalar gasM_;       

    //- Old time gas mass [kg]
    scalar gasM0_; 
    
    //- Species mass fractions [-]
    scalarList Y_;  
    
    //- Old time mass fractions [-]
    scalarList Y0_;
    
    //- molarFraction [-]
    scalarList M_;
    
    //- Old time molarFraction [-]
    scalarList M0_;        

    //- Correct the mass fractions to sum to 1
    virtual void correctMassFractions();
    
    //- Correct the molar fractions to sum to 1
    virtual void correctMolFractions();

public:           

    //- Runtime type information
    TypeName("none");    
    
    // Declare run-time constructor selection table

    declareRunTimeSelectionTable
    (
        autoPtr,
        gapGasModel,
        dictionary,
        (

            const fvMesh& mesh,
            const materials& mat,
            const volScalarField& T,
            const volVectorField& DorDD,
            const fissionGasRelease& fgr,
            const dictionary& gapGasDict
        ),
        (mesh, mat, T, DorDD, fgr, gapGasDict)
    );

    // Constructors

    //- Construct from mesh, materials and dict
    gapGasModel
    (    
        const fvMesh& mesh,
        const materials& mat,
        const volScalarField& T,
        const volVectorField& DorDD,
        const fissionGasRelease& fgr,
        const dictionary& gapGasDict
    );

    //- Return a pointer to a new thermoMechanicalProperties
    static autoPtr<gapGasModel> New
    (   
        const fvMesh& mesh,
        const materials& mat,
        const volScalarField& T,
        const volVectorField& DorDD,
        const fissionGasRelease& fgr,
        const dictionary& solverDict
    );

    
    //- Destructor
    virtual ~gapGasModel();


    // Access Functions         
                    
    //- Return gas pressure
    inline scalar p() const
    {
        return gasP_;
    }  
    
    //- Return gap temperature
    inline scalar T() const
    {
        return gapT_;
    }  
    
    //- Return gas mass
    inline scalar m() const
    {
        return gasM_;
    }    
    
    //- Return old time gas mass
    inline scalar m0() const
    {
        return gasM0_;
    }          
    
    //- Return mass fractions
    inline scalarList Y() const
    {
        return Y_;
    }        
    
    //- Return old time mass fractions
    inline scalarList Y0() const
    {
        return Y0_;
    }           
    
    //- Return molar fractions
    inline scalarList M() const
    {
        return M_;
    }        

    //- Return old time molar fractions
    inline scalarList M0() const
    {
        return M0_;
    }          
    
    //- Return gas component names
    inline wordList species() const
    {
        return species_;
    }          
    
    //- Return gas component molecular weights
    inline scalarList speciesW() const
    {
        return speciesW_;
    }            
    
    // Member Functions

    //- Gap gas thermal conductivity
    virtual scalar kappa(const scalar T) const
    {
        FatalErrorIn("Foam::gapGasModel::kappa(const scalar T)")
        << "function called for \"none\" gas model. Select a gap model." 
        << endl << Foam::abort(FatalError);

        return 0;
    };
    
    //- Return the accommodation coefficient
    virtual scalar a(const scalar T) const 
    {
        FatalErrorIn("Foam::gapGasModel::a(const scalar T)")
        << "function called for \"none\" gas model. Select a gap model." 
        << endl << Foam::abort(FatalError);

        return 0;
    };
    
    //- Modify mass and mass fractions based on fgr
    virtual void correctMass();

    //- Correct the gas thermodynamic properties
    virtual void correct() {}; 

    //- Update variables at the end of time step
    virtual void updateVariables(); 
        
    // IO
    
    //- Read data from time directory
    virtual bool readData(Istream& is)
    {
            return !is.bad();
    }
    
    //- Write data to time directory
    virtual bool writeData(Ostream& os) const 
    {
            return os.good();
    };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
