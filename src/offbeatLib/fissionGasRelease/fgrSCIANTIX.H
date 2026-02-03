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
    Foam::fgrSCIANTIX

Description
	This class integrates the 0-D SCIANTIX code (developed at Polimi) into OFFBEAT 
	as the fission gas release (FGR) model. The coupling process between OFFBEAT and 
	SCIANTIX proceeds as follows:

	1. At each outer iteration and for every cell of fuel-type materials (such as `UO2`), 
	OFFBEAT passes the necessary data to SCIANTIX for calculations. This data set includes:
		- History variables, i.e., local oldTime and current stresses, temperature, and 
		fission rate.
		- Previous time step values of the SCIANTIX variables for the specific cell.

	2. Using the previous time step values as initial conditions, SCIANTIX computes 
	the updated values for its internal variables.

	3. At the end of each outer iterations, certain updated values calculated by 
	SCIANTIX, such as the FGR fraction or gaseous swelling, are passed to OFFBEAT.

	4. At the end of each time step, steps 1, 2, and 3 are repeated. This final pass 
	stores the updated local SCIANTIX variables into the corresponding OFFBEAT fields, 
	completing the synchronization process.
	
________________________________________________________________________________

\par Options	
	The `fgrSCIANTIX` model requires the user to specify a few additional parameters in 
	the `fgrOptions` sub-dictionary located inside the `solverDict`. Additionally,
    SCIANTIX requires the user to set the specific SCIANTIX input settings. See the next 
    section for more information or consult SCIANTIX documentation.

    \warning
    There are two ways to provide SCIANTIX input settings when using OFFBEAT. 
    One option is to define these settings in the `input_settings.txt` file within 
    the OFFBEAT case folder (i.e., at the same level as the `constant` or 
    `system` folders). However, note that if OFFBEAT is launched from outside the 
    case folder using the `-case` argument, SCIANTIX assumes the `input_settings.txt` 
    file path is the same as the folder from which the execution command is launched.
    <br><br>
    To address this, a new input reading option has been added to OFFBEAT. SCIANTIX 
    input settings can now be defined in a sub-dictionary of `fgrOptions` within `solverDict`. 
    If the code detects this dictionary, it reads the input from there; otherwise, 
    it defaults to the original approach, reading from the `input_settings.txt` file. 
    In both cases, an `input_check.txt` file is created for input validation. An example 
    of the new input format can be found in `Cases/testCases/generic_pwr_1D_newSCIANTIXinput/`.
    
    Parameters in `fgrOptions`:
	
	- <b>`relax`</b> - Relaxation factor applied to inter/intra-gaseous swelling and fgr
	fraction. It can improve convergence when dealing with strongly coupled
	simulations. **It is set to 1 by default**.
	- <b>`nFrequency`</b> - An integer defining how many outer iteration to skip before
	updating the OFFBEAT/SCIANTIX coupling. It can speed up the simulation (i.e.
	less SCIANTIX runs) but it should not be too large, as one risks not attaining
	proper convergence. **It is set to 1 by default**.
	- <b>`addToRegistry`</b> - Keyword to switch on/off the creation of the SCIANTIX fields
	as separate volumetric fields. **It is set to `off` by default**.
	
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
	coupling between OFFBEAT and SCIANTIX .
	
	\verbatim
		fgr     SCIANTIX;

		fgrOptions
		{
			//- relax is 1 by default
			relax  1;

			//- nFrequency is 1 by default
			nFrequency  1;

			//- addToRegistry is off by default
			addToRegistry  off;
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
	
	Additionally, the code SCIANTIX requires its own set of inputs to be
	properly set. One way to set the SCIANTIX parameters is to properly define in the
    case folder the `input_settings.txt` file (see SCIANTIX documentation for more details):
	
	\verbatim
		0 #  verification (0 = no verification)
		1 #  grain growth (0 = none; 1 = ainscough)
		1 #  inert gas behavior (0 = none; 1 = do IGB)
		1 #  gas diff coeff (0 = 0 coeff (trial); 1 = Turnbull, 2 = Matzke (with small
		OFFBEAT custom addition) )
		1 #  intra bbl evo (1 = Pizzocri et al., 2018; 2 = White, Tucker 1983)
		1 #  intra bubble_radius (1 = Olander&Wongy, 2006, 2= Ronchi 1972; 3 = Dollings
		1977; 4 = Losonen 2002; 5 = Spino 2005; )
		1 #  re-solution (1 = Turnbull 1971; 2 = Oldanter 1976)
		1 #  trapping (1 = Ham 1958)
		1 #  nucleation (1 = Baker 1971)
		1 #  DiffSolver (1 = SDA, Pizzocri, 2 = FORMAS)
		1 #  format_out (1 = output.txt, values separated by tabs)
		1 #  gb vac diff coeff (1 = Reynolds and Burton, 1979; 2 = Pastore 2015)
		1 #  gb behavior (1 = Pastore et al., 2013; Barani et al., 2017)
		0 #  gb micro-cracking (0 = none; 1 = Barani et al., 2017)
		0 #  hbs formation (0 = non active)
		0 #  fuel/reactor couple for burnup calculations (0 = UO2/PWR)
		1 #  gas effective (=0) or single atom (=1) diffusion coefficient - OFFBEAT
		custom addition
		0 #  add boundary gas sweeping (1=on) - OFFBEAT custom addition
		0 #  micro cracking model span parameter (0 = default value of 10; 1 = original
		value SCIANTIX of 10; 2 = Value from MOOSE of 5)
	\endverbatim

    An alternative (which take priority in OFFBEAT) is to define the SCIANTIX
    options directly in the `fgrOptions` dictionary.

    \verbatim
        fgrOptions
        {
            nFrequency  1;
            relax       1;

            SCIANTIX
            {
                iverification                                   0; // (0= no verification)
                igrain_growth                                   1; // (1 = ainscough)
                iinert_gas_behavior                             1; // (1= do IGB)
                igas_diffusion_coefficient                      1; // (1= Turnbull et al., 1988)
                iintra_bubble_evolution                         1; // (1=Pizzocri et al., 2018)
                ibubble_radius                                  1; // (1= Olander&Wongy, 2006)
                iresolution_rate                                1; // (1=Turnbull 1971)
                itrapping_rate                                  1; // (1= Ham 1958)
                inucleation_rate                                1; // (1= Baker 1971)
                isolver                                         1; // (1= SDA, Pizzocri et al., 2019)
                iformat_output                                  1; // (1 = output.txt, values separated by tabs)
                igrain_boundary_vacancy_diffusion_coefficient   1; // (1= Reynolds and Burton, 1979)
                igrain_boundary_behaviour                       1; // (1= do InterGranularGasBehavior - Pastore et al., 2013; Barani et al., 2017)
                igrain_boundary_micro_cracking                  1; // (1 = Barani et al., 2017)
                igrain_recrystallization                        0; // (0 = non active)
                ifuel_reactor_type                              0; // (0=UO2/PWR)
                igas_effective_coefficient                      1; // gas effective (=0) or single atom (=1)
                igas_sweeping                                   1; // add boundary gas sweeping
            }

        }
	\endverbatim
    
SourceFiles
    fgrSCIANTIX.C

\todo
    - use cutting plane or fixed 1D geometry to solve SCIANTIX only at few
    axial locations?  

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)\n
    A. Gianesello - INSA, Lyon

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef fgrSCIANTIX_H
#define fgrSCIANTIX_H

#include "fissionGasRelease.H"
#include "volFields.H"
#include "List.H"
#include "globalOptions.H"

#include "MainSCIANTIX.h"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class fgrSCIANTIX Declaration
\*---------------------------------------------------------------------------*/

class fgrSCIANTIX
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

        //- Rod is refabricated, so fgr is a % of start of transient inventory
        bool refabricatedRod_;

        //- Dictionary for SCIANTIX input
        dictionary inputSCIANTIX_;
        
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

        scalarField& Intragranular_bubble_concentration_;
        scalarField& Intragranular_bubble_radius_;
        scalarField& Intergranular_bubble_concentration_;
        scalarField& Intergranular_atoms_per_bubble_;
        scalarField& Intergranular_vacancies_per_bubble_;
        scalarField& Intergranular_bubble_radius_;
        scalarField& Intergranular_bubble_area_;
        scalarField& Intergranular_bubble_volume_;
        scalarField& Intergranular_fractional_coverage_;
        scalarField& Intergranular_saturation_fractional_coverage_;
        scalarField& Intergranular_fractional_intactness_;
        scalarField& timeStepNumber_;

        PtrList<scalarField> grainModes_;

        //- Value of current time
        scalar time_;        

        //- Value of current outerIteration
        int outerIteration_;        

        // How many outerIterations should be completed for SCIANTIX to run
        int nFrequency_;      

        //Relaxation factor
        scalar relax_;           

        //- Pointer to fields often used by the class
        const volScalarField* Q_;
        const volScalarField* Bu_;
        const volScalarField* T_;
        const volScalarField* rho_;
        const volSymmTensorField* sigma_;

        //- True if input has been read
        bool inputRead_;

        //- % FGR (used for functionObject)
        scalar fgrPercent_;

        //- FGR (m3) in the 360 deg rod @ T= 293 K and P=101325 Pa
        scalar fgrM3_;

        //- FG produced (m3) in the 360 deg rod @ T= 293 K and P=101325 Pa
        scalar fgpM3_;

    // Private Member Functions

        //- Read input sciantix
        virtual void readSCIANTIXinput();

        //- Disallow default bitwise copy construct
        fgrSCIANTIX(const fgrSCIANTIX&);

        //- Disallow default bitwise assignment
        void operator=(const fgrSCIANTIX&);

public:

    //- Runtime type information
    TypeName("SCIANTIX");


    // Constructors

        //- Construct from mesh, materials and dict
        fgrSCIANTIX
        (
          const fvMesh& mesh,
          const materials& mat,
          const dictionary& fissionGasReleaseDict
        );


    //- Destructor
    ~fgrSCIANTIX();


    // Acess Functions

        //- Returns %fgr (for functionObject)
        virtual const scalar fgrPercent() const
        {
            return fgrPercent_;
        }

        //- Returns fgr (m3) (for functionObject) @ T= 293 K and P=101325 Pa
        virtual const scalar fgrM3() const
        {
            return fgrM3_;
        }

        //- Returns fg produced (m3) (for functionObject) @ T= 293 K and P=101325 Pa
        virtual const scalar fgpM3() const
        {
            return fgpM3_;
        }

    // Member Functions
        
       //- Update fission gas distribution in fuel grains
       virtual void correct();

       //- Pass previous time-step cell variables to SCIANTIX
       void setSCIANTIXVariables(const label);

       //- Pass new valus of variables to OFFBEAT and updates volFields
       void updateSCIANTIXVariables(const label);
        
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