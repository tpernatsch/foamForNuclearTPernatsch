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
    Foam::fissionGasRelease

Description
	As the name implies, the `none` fission gas release class in OFFBEAT allows you
	to switch off the fgr models and have a simulation where fgr production, release
	and related gaseous swelling are not considered. This means that no additional 
  fields are created.
	
	\warning
	Although with the `none` fgr class no additional fields are created for the
	fission gases itself (e.g. concentration, release fraction etc.) the basic fgr
	class generates the following fields:	<b>`grainRadius`</b>, <b>`oxygenMetalRatio`</b>
	as well as <b>`interGranularSwelling`</b> and <b>`interGranularSwelling`</b>.
  While this specific aspect might change in future version of OFFBEAT, at the
	moment selecting `none` means that the previously mentioned fields will be
	created and will remain constant during the simulation.
	
________________________________________________________________________________

Usage
	Here is an example of (part of) the `solverDict` file to be used to activate the
	`none` fissionGasRelease class
	
	\verbatim
    fgr none;
	\endverbatim

SourceFiles
    fissionGasRelease.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)\n

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef fissionGasRelease_H
#define fissionGasRelease_H

#include "primitiveFields.H"
#include "volFields.H"
#include "IOdictionary.H"
#include "materials.H"
#include "fuelMaterial.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class fissionGasRelease Declaration
\*---------------------------------------------------------------------------*/

class fissionGasRelease
{
    // Private data
  

    // Private Member Functions

        //- Disallow default bitwise copy construct
        fissionGasRelease(const fissionGasRelease&);

        //- Disallow default bitwise assignment
        void operator=(const fissionGasRelease&);

protected:
    
    // Protected data

        //- Reference to object fvMesh
        const fvMesh& mesh_;

        //- Const reference to the materials class
        const materials& mat_;

        //- Fgr options dictionary
        const dictionary fissionGasReleaseDict_;

        //- Current fgr fraction (in percent). Necessary for time-stepping
        scalarField fgr_;

        //- Old fgr fraction (in percent). Necessary for time-stepping
        scalarField fgr0_; 

        //- Additional mole of Xe and Kr released during current time step
        scalar molFGR_;

        //- Additional mole of He released during current time step
        scalar molHe_;   

        //- Maximum allowed local change in fgr (%) during time-step 
        //  (for setting next deltaT value)
        scalar maxLocalDeltaFgr_;

        //- Maximum allowed local change in fgr (%) during time-step 
        //  (for setting next deltaT value)
        scalar maxTotalDeltaFgr_;    

        //- Maximum allowed number of moles released in the model during time-step 
        //  (for setting next deltaT value) NOTE: deprecated method
        scalar maxFGR_;

        //- TimeStepping criterion (either use maxFGR_ or use maxTotal/LocalFgrChange)
        //  NOTE: This is only temporary, it will be removed once the maxFGR_ method is removed
        word timeSteppingCriterion_;

        //- Grain radius
        volScalarField& grainRadius_;

        //- OxygenMetal (O/M) ratio
        volScalarField& oxygenMetalRatio_;

        //- Intra-granular gas swelling        
        volSymmTensorField intragranularGasSwelling_;
        
        //- Inter-granular gas swelling        
        volSymmTensorField intergranularGasSwelling_;

public:

    //- Runtime type information
    TypeName("none");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            fissionGasRelease,
            dictionary,
            (
              const fvMesh& mesh,
              const materials& mat,
              const dictionary& fissionGasReleaseDict
            ),
            (mesh, mat, fissionGasReleaseDict)
        );

    // Constructors

        //- Construct from mesh, materials and dict
        fissionGasRelease
        (
          const fvMesh& mesh,
          const materials& mat,
          const dictionary& fissionGasReleaseDict
        );


    // Selectors

        //- Select
        static autoPtr<fissionGasRelease> New
        (
          const fvMesh& mesh,
          const materials& mat,
          const dictionary& solverDict
        );


    //- Destructor
        virtual ~fissionGasRelease();


    // Access Functions
     
        const scalar& molFGR() const
        {
          return molFGR_;
        };

        const scalar& molHe() const
        {
          return molHe_;
        };


    // Member Functions

        virtual void correct() {};
          
        //- Update variables at the end of time step
        virtual void updateVariables(){}; 
         
        //- Return ordered list of released gas mixture component names. 
        //- Necessary for supplying the mass source to the gap gas model
        virtual const wordList gasComponents() const
        {
          const wordList emptyList;

          return emptyList;
        };
         
        //- Return ordered list of released gas mixture component mols.
        //- The mols are the additional mols released during current time step.
        //- Necessary for supplying the mass source to the gap gas model
        virtual const scalarList gasMols() const
        {
          const scalarList emptyList;

          return emptyList;
        };

        virtual scalar nextDeltaT(); 
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
