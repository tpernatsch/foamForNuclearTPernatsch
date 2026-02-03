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
    Foam::fgrSCIANTIXRIA

Description
	This class is designed for use in RIA simulations, specifically after a restart.
    The intended workflow is for the user to apply the 'SCIANTIX' fission gas release 
    model during the base irradiation phase and then switch to 'SCIANTIXRIA' for the 
    restart or RIA simulation.
    
    In this class, the SCIANTIX module itself is not actively called. Instead, the 
    SCIANTIX fields - representing grain and grain boundary gas concentrations - 
    computed during the base irradiation are read from the initial restart time folder.

    During each iteration, the 'SCIANTIXRIA' class checks whether a given fuel cell has 
    exceeded certain conditions:
    
    - The burnup or temperature threshold for HBS, if HBS release is activated. In this case
    the grain and grain boundary gases are vented (added to the release fraction) 
    and the intra- and inter-granular swelling fields are set to 0.
    - A damage threshold where the 'damage' field is typically calculated using the 
    'damageModel' class in the constitutive law behavior. In this case, only the 
    inter-granular gases are release and the corresponding swelling is set to 0
	
________________________________________________________________________________

\par Options	
	The `fgrSCIANTIXRIA` model requires the user to specify a few additional parameters in 
	the `fgrOptions` sub-dictionary located inside the `solverDict`. 
    
    Parameters in `fgrOptions`:
	
	- <b>`relax`</b> - Relaxation factor applied to inter/intra-gaseous swelling and fgr
	fraction. It can improve convergence when dealing with strongly coupled
	simulations. **It is set to 1 by default**.
	- <b>`releaseHBS`</b> - Switch to activate HBS release model. **It is set to 
    true by default**.
	- <b>`buReleaseThresholdHBS`</b> - Burnup HBS release threshold in MWd/t. 
    **It is set to 80000 by default**.
	- <b>`temperatureReleaseThresholdHBS`</b> - Temperature HBS release threshold 
    in K. **It is set to 1000 by default**.
	- <b>`damageReleaseThreshold`</b> - Damage release threshold in fraction.
    **It is set to 0.85 by default**.
	
	\note
	To minimize the creation of multiple volumetric fields (such as `volScalarFields` 
	or `volVectorFields`), SCIANTIX fields are, by default, defined as standard fields 
	(`scalarFields` or `vectorFields`). This is achieved using a template function located 
	in `globalFields.H`.
	<br><br>
	Typically, these fields are read from a `regIOObject` named `SCIANTIXfields`, which 
	is saved at each write time step. At the top of this dictionary, a summary of the 
	current fission gas release (FGR) is included. If a specific field is missing from the 
	dictionary, it is initialized with a default value. However, if an `IOObject` 
	corresponding to one of these fields exists in the initial time step folder, the 
	field is read from that file. This approach can be helpful when performing advanced 
	operations, such as using OpenFOAM’s mapping routines.
	<br><br>
	In certain scenarios, adding fields to the registry is necessary for specific tasks, 
	such as visualization in ParaFoam, probing, or mapping. This behavior is controlled by 
	the <b>`addToRegistry`</b> keyword in the FGR options dictionary, enabling users to specify 
	whether fields should be registered as volumetric fields.
	
________________________________________________________________________________

Usage
	Here is an example of (part of) the `solverDict` file to be used to activate the
	SCIANTIXRIA class.
	
	\verbatim
		fgr     SCIANTIXRIA;

		fgrOptions
		{
			//- relax is 1 by default
			relax  1;

            // Switch to activate HBS release model
            releaseHBS true;

            // Bu HBS release threshold in MWd/t
            buReleaseThresholdHBS 80000;

            // Temperature HBS release threshold in K
            temperatureReleaseThresholdHBS 1000;

            // Damage release threshold in fraction
            damageReleaseThreshold 0.85;
		}
	
		//... other physics subdictionaries

		materials
		{
			fuel
			{
				type UO2;

				//... other material specific options

				rGrain 5e-6;

			}
			//... remaining material subdictionaries
		}
	\endverbatim
	
	There is no need to set the SCIANTIX input file for this class as SCIANTIX
    code is not actually called.
    
SourceFiles
    fgrSCIANTIXRIA.C

\todo
    - use cutting plane or fixed 1D geometry to solve SCIANTIX only at few
    axial locations?  

\mainauthor
    M. Reymond - EPFL/CEA
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Vanyi - EPFL

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef fgrSCIANTIXRIA_H
#define fgrSCIANTIXRIA_H

#include "fissionGasRelease.H"
#include "volFields.H"
#include "List.H"
#include "globalOptions.H"

#include "MainSCIANTIX.h"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class fgrSCIANTIXRIA Declaration
\*---------------------------------------------------------------------------*/

class fgrSCIANTIXRIA
:
    public fissionGasRelease,
    public regIOobject
{
    // Private data
        
        //- List of zone names where to apply the model
        List<word> zoneNames_;

        //- regIOObject for reading/writing SCIANTIX fields
        dictionary sciantixFieldsDict_;

        //- Add sciantix fields to registry (write to disc)
        bool addToRegistry_;
        
        //- Additional dependencies needed by the SCIANTIX code
        scalarField& Gas_produced_;
        scalarField& Gas_grain_;
        scalarField& Gas_grain_solution_;
        scalarField& Gas_grain_bubbles_;
        scalarField& Gas_boundary_;
        scalarField& Gas_released_;
        scalarField& Gas_released_t0_;
        scalarField& Effective_burn_up_;
        scalarField& Helium_produced_;
        scalarField& Helium_grain_;
        scalarField& Helium_grain_solution_;
        scalarField& Helium_grain_bubbles_;
        scalarField& Helium_boundary_;
        scalarField& Helium_released_;

        //- Relaxation factor
        scalar relax_;

        //- True if input has been read
        bool inputRead_;

        //- Activate release due to HBS formation
        Switch releaseHBS_;

        //- In case releaseHBS_ is active, then select temperature threshold
        // for release
        scalar buReleaseThresholdHBS_;

        //- In case releaseHBS_ is active, then select burnup threshold
        // for release
        scalar temperatureReleaseThresholdHBS_;

        //- Damage threshold for overfragmention release
        scalar damageReleaseThreshold_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        fgrSCIANTIXRIA(const fgrSCIANTIXRIA&);

        //- Disallow default bitwise assignment
        void operator=(const fgrSCIANTIXRIA&);

public:

    //- Runtime type information
    TypeName("SCIANTIXRIA");


    // Constructors

        //- Construct from mesh, materials and dict
        fgrSCIANTIXRIA
        (
          const fvMesh& mesh,
          const materials& mat,
          const dictionary& fissionGasReleaseDict
        );


    //- Destructor
    ~fgrSCIANTIXRIA();


    // Member Functions
        
       //- Update fission gas distribution in fuel grains
       virtual void correct();
        
       //- Update variables at the end of time step
       virtual void updateVariables();
       
       //- Return ordered list of released gas mixture component names. 
       //- Necessary for supplying the mass source to the gap gas model
       virtual const wordList gasComponents() const;
       
       //- Return ordered list of released gas mixture component mols.
       //  The mols are the total released mols since the start of the simulation.
       //  Necessary for supplying the mass source to the gap gas model
       virtual const scalarList gasMols() const;            
           

    // IO
        //- Write data to time directory
        virtual bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif