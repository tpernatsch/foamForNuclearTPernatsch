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
    Foam::gapTRISO

Description
	Model for gap gas composition, volume, temperature and pressure for TRISO.

    #### Formulation

	Similar to the `Frapcon` model in terms of main assumptions and algorithm used 
	to calculate the gap volume. Check `Frapcon` documentation for more 
	information.
	
	\warning 
	The volumes calculated and printed are in m3, and refer to the modeled geometry.
	
	In addition the gapTRISO model calculates the production of CO using
	either the GA or Proksch models.
	
	For the GA model, one has that:
	
	$$
	\begin{aligned}
	CO_{prod}=min(0.625, 1.64e^{\frac{-3311}{T}})
	\end{aligned}
	$$
	
	where:

	- \f$CO_{prod}\f$ is produced CO per fission(atoms/fission),
	- \f$T\f$ is temperature(K),
	- 0.625 is the limit set by Parfume.
	
	For the Proksch model, one has that:
	
	$$
	\begin{aligned}
	log\frac{O/F}{t^2}=-0.21-\frac{8500}{T}
	\end{aligned}
	$$
	
	where:
	
	- \f$O/F\f$ is oxygen release at the end of irradiation(atoms/fission),
	- \f$t\f$ is irradiation time(days),
	- \f$T\f$ is time-averaged particle surface temperature(K).
	
	The maximum value of \f$O/F\f$ is limited by the following expression:
	
	$$
	\begin{aligned}
	(O/F)_{max}=0.40f_U+0.85f_{Pu}
	\end{aligned}
	$$
	
	where \f$f_U\f$ and \f$f_{Pu}\f$ are the fission fractions of \f$^{235}U\f$ and
	\f$^{239}Pu + ^{241}Pu\f$, respectively.
	
	We assume that all the oxygen forms carbon monoxide, then \f$CO_{prod}=O/F\f$.
	
________________________________________________________________________________

\par Options
    The `TRISO` class in `OFFBEAT` requires the user to set various parameters
    in `gapGasOptions` (a subdictionary of the `solverDict`) and `gapGas` 
    (a dictionary in the `0/uniform/` folder), as described below.

    Parameters in `gapGasOptions`:

    - <b>`heatSourceName`</b> - Name of heat source field. **It is `Q` by default**.
    - <b>`COProductionModel`</b> - CO production model, either `Proksch` (default) or `GA`.
    - <b>`model3D`</b> - Is the model in 3D? Set to `false` by default.
    - <b>`bufferPorosity`</b> - The porosity of `Buffer`.
    - <b>`bufferZones`</b> - Name of the zones corresponding to `Buffer` material.
    - <b>`gapPatches`</b> - Specifies the patches for the fuel-clad gap (e.g., 
    `fuelOuter cladInner`).
    - <b>`gapVolumeOffset`</b> - Offset for the initial gap volume (by default: `0.0`).
    - <b>`gasReserveVolume`</b> - Initial volume of reserve gas for gap pressure 
    calculation (by default: `0.0`).
    - <b>`gasReserveTemperature`</b> - Initial temperature of reserve gas, in Kelvin 
    (by default: `290`). Alternatively the user can provide a time-dependent list 
    using the following keyword
    - <b>`gasReserveTemperatureList`</b> - If present it is read in place of 
    `gasReserveTemperature`. It is a sub-dictionary used to provide a time-dependent table 
    of reserve temperatures.
    - <b>`type`</b> - Parameter of `gasReserveTemperatureList` subdictionary; it specifies the 
    table type (`table`).
    - <b>`values`</b> - Parameter of `gasReserveTemperatureList` subdictionary; it is a list of 
    time-value pairs for temperature (e.g., `[(0.0, 290), (100.0, 590)]`).
    - <b>`outOfBounds`</b> - Parameter of `gasReserveTemperatureList` subdictionary; it handles 
    the method for out-of-bounds values 
    (default: `clamp`).
    - <b>`interpolationScheme`</b> - Parameter of `gasReserveTemperatureList` subdictionary; 
    it selects the interpolation method between time points (e.g., `linear`).
    - <b>`includeCentralHole`</b> - Toggles inclusion of central hole in the fuel 
    (default: `on`).
    - <b>`includeDishes`</b> - Toggles inclusion of dish geometry (default: `on`).

    Parameters in `gapGas`:

    - <b>`massFractions`</b> - Specifies the initial mass fractions of gases 
    (`Ar`, `He`, `Kr`, `Ne`, `Rn`, `Xe`, `CO`).
    - <b>`gasPressureType`</b> - Defines the method for gas pressure calculation 
    (options: `fromModel`, `fixed`, or `fromList`).
    - <b>`gasPressure`</b> - Fixed gas pressure value when `gasPressureType` is set 
    to `fixed`.
    - <b>`gapPressureList`</b> - Activated only if `gasPressureType` is set as 
    `fromList`, this is a sub-dictionary used to provide a time-dependent table 
    for gas pressure.
    - <b>`type`</b> - Parameter of `gapPressureList` subdictionary; it specifies the 
    table type (`table`).
    - <b>`values`</b> - Parameter of `gapPressureList` subdictionary; it is a list of 
    time-pressure pairs (e.g., `[(0.0, 1e5), (100.0, 1e6)]`).
    - <b>`outOfBounds`</b> - Parameter of `gapPressureList` subdictionary; it handles the
    method for out-of-bounds values (default: `clamp`).
    - <b>`interpolationScheme`</b> - Parameter of `gapPressureList` subdictionary; it specifies the 
    interpolation method between time points (e.g., `linear`).

________________________________________________________________________________

Usage	
	Here is a code snippet of the `solverDict` to be used for activating the
	`Frapcon` gapGas model:
	
	\verbatim
		gapGas                  TRISO;

		//... rest of solverDict
		
		// Option subdict for gapGas model
		gapGasOptions
		{
			gapPatches ( bufferOuter ipycInner );
			bufferOuterPatches bufferOuter;
			bufferPorosity  0.4054;
			bufferZones     Buffer;
			COProductionModel  Proksch;
			model3D    false;
			gapVolumeOffset 0.0;
			gasReserveVolume 0.0;
			gasReserveTemperature 290;

			// // Alternatively one can provide a time dependent list
			// gasReserveTemperatureList
			// {
			//      type    table;
			//
			//        values
			//        (
			//            (0.0    290)
			//            (100.0  590)
			//        );
			//
			//        outOfBounds         clamp;     // optional out-of-bounds handling
			//        interpolationScheme linear;    // optional interpolation method
			// }
		}
	\endverbatim
	
	In `gapGas file` inside the time step folder (e.g. 0/uniform or 0/):
	
	\verbatim
		massFractions
		{
			Ar 0;
			He 0;
			Kr 0;
			Ne 0;
			Rn 0;
			Xe 0;
			CO 0;
		};

		gasPressureType                 fromModel;// fixed;// fromList;
		gasPressure                     0;

		//  gapPressureList
		//  {
		//  type    table;
		//
		//      values
		//      (
		//          (0.0    1e5)
		//          (100.0  1e6)
		//      );

		//      outOfBounds         clamp;     // optional out-of-bounds handling
		//      interpolationScheme linear;    // optional interpolation method
		//  }
	\endverbatim

SourceFiles
    gapTRISO.C

\mainauthor
    F. Xiang - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
    Laboratory for Reactor Physics and Systems Behaviour) & XJTU (Xi'an JiaoTong
    University, China, Nuclear Thermal-hydraulic Laboratory)

\contribution
    A. Scolaro - EPFL
    E. Brunetto - EPFL

\date
    Feburary 2023

\*---------------------------------------------------------------------------*/

#ifndef gapTRISO_H
#define gapTRISO_H

#include "gapGasModel.H"
#include "Function1.H"
#include "globalOptions.H"
#include <set>
#include "Table.H"
#include "polyTopoChanger.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class gapTRISO Declaration
\*---------------------------------------------------------------------------*/

class gapTRISO
:
    public gapGasModel
{
    // Private data

    // Private Member Functions

        //- Calculate the initial gas mass
        virtual void calcInitialMass();

        //- Return the mixture density calculated with ideal gas law
        virtual scalar rhoMixture(const scalar&, const scalar&) const;

        //- Find patch IDs for patches grouped under given name
        virtual void findPatchIDs(const word groupName);

        //- Estimate the current gap quantities (volume, temperature, VoverT)
        virtual void correctGap();

        //- Estimate the volume of gas inside the Buffer
        virtual void correctBuffer();

        //- Override the same functions in gapGasModel.
        virtual void correctMassFractions();
        virtual void correctMolFractions();
        virtual void correctMass();

        // //- Update variables at the end of time step
        // virtual void updateVariables();

        //- To calculate the production of CO in current time step (moles)
        virtual void productionCO();

        //- production rate of CO, GA model (atoms/fission)
        virtual scalar GAmodel();

        //- production rate of CO, Proksch model (atoms/fission)
        virtual scalar ProkschModel();


        //- Disallow default bitwise copy construct
        gapTRISO(const gapTRISO&);

        //- Disallow default bitwise assignment
        void operator=(const gapTRISO&);

protected:

    // Protected data
        const word heatSourceName_;

        const word CO_method_;
        const bool is3D_;

        //- Cracks volume [m3]
        scalar crackV_;

        //- Gap volume offset [m3]
        scalar gapOffsetV_;

        //- Additional gas reserve volume [m3]
        scalar reserveV_;

        //- The porosity of Buffer
        scalar porosity_;

        //- Volume of gas inside the Buffer [m3]
        scalar bufferV_;

        //- Crack gas average temperature [K]
        scalar crackT_;

        //- Additional gas reserve temperature [K]
        scalar reserveT_;

        //- Temperature of gas inside the Buffer [m3]
        scalar bufferT_;

#ifdef OPENFOAMFOUNDATION
        //- Additional gas reserve temperature (time list) [K]
        autoPtr<Function1s::Table<scalar>> reserveTlist_;  
#elif OPENFOAMESI        
        //- Additional gas reserve temperature (time list) [K]
        autoPtr<Function1<scalar>> reserveTlist_;   
#endif 

        //- Gap sum(Vi/Ti) for each face of the gap [m3/K]
        scalar gapVoverT_;

        //- Crack sum(Vi/Ti) for each fuel cell [m3/K]
        scalar crackVoverT_;

        //- Buffer total_V/average_T [m3/K]
        scalar bufferVoverT_;

        //- Value of current time
        scalar currentTime1_;

        //- Value of current time
        scalar currentTime2_;

        //- Value of total fission (fission)
        scalar total_fission_;

        //- read the total fission (fission) of previous simulation (for restart)
        scalar total_fission0_;

        scalar Tint_;

        //- read the time integrated temperature of previous simulation (for restart)
        scalar Tint0_;

        //- Increment of CO [moles]
        scalar Mco_;

        //- List of scaling factors (one list per patch)
        List<scalarField> scalingFactors_;

        //- Hash table of patch groups bounding the gap (gap, fuelTop, fuel bottom)
        //  TODO: separate fuelGap and claddingGap
        autoPtr<Foam::HashTable<std::set<label>, Foam::word>> patchIDsPtr_;

        //- Read As parameters for Frapcon conductivity model
        scalarList conductivityAs_;

        //- Read Bs parameters for Frapcon conductivity model
        scalarList conductivityBs_;

        //- Type of pressureModel: calculated or read from a list?
        word gasPType_;

#ifdef OPENFOAMFOUNDATION
        //- Time dependent list of gas pressure values
        autoPtr<Function1s::Table<scalar>> gasPressureList_;   
#elif OPENFOAMESI        
        //- Time dependent list of gas pressure values
        autoPtr<Function1<scalar>> gasPressureList_;   
#endif  

public:

    //- Runtime type information
    TypeName("TRISO");


    // Constructors

        //- Construct from mesh, materials and dict
        gapTRISO
        (
            const fvMesh& mesh,
            const materials& mat,
            const volScalarField& T,
            const volVectorField& DorDD,
            const fissionGasRelease& fgr,
            const dictionary& gapGasDict
        );


    //- Destructor
    ~gapTRISO();


    // Member Functions

        //- Gap gas thermal conductivity
        //- (CO, Ne from Matpro, other gas from Frapcon manual)
        virtual scalar kappa(const scalar T) const;

        //- Return the accommodation coefficient (from Frapcon manual)
        virtual scalar a(const scalar T) const;

        //- Correct the gas thermodynamic properties
        virtual void correct();


    // IO

        //- Read data from time directory
        virtual bool readData(Istream& is);

        //- Write data to time directory
        virtual bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
