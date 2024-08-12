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
    Model for gap gas composition, volume, temperature and pressure for TRISO

    \htmlinclude gapGasModel_gapTRISO.html

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

        //- Estimate the volume of gas inside the buffer
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

        //- The porosity of buffer
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

        //- Read As parameters for FRAPCON conductivity model
        scalarList conductivityAs_;

        //- Read Bs parameters for FRAPCON conductivity model
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
            const dictionary& gapGasOptDict
        );


    //- Destructor
    ~gapTRISO();


    // Member Functions

        //- Gap gas thermal conductivity
        //- (CO, Ne from MATPRO, other gas from FRAPCON manual)
        virtual scalar kappa(const scalar T) const;

        //- Return the accommodation coefficient (from FRAPCON manual)
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
