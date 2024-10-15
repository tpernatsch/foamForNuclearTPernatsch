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
    Foam::gapFRAPCON

Description
    Model for gap gas composition, volume, temperature and pressure. 
    Some routines are derived from FRAPCON manual. 

    NOTE: the volumes are in m3, and refer to the modeled geometry.
    e.g. consider a full rod plenum volume of 36 cm3. If the rod is reproduced 
    with a 2-degree wedge mesh, the plenum volume in the model is:
    \f[
                    V_model_plenum = 36/1000/1000/360*2 m3 = 2E-07 m3
    \f]  

    The main assumptions are:
    - ideal gas law
    - each region of the free volume has different V and T
    - the pressure is instantaneously equalized everywhere in the free volume
    
    The gas pressure is going to be calculated as:
    \f[
                    p = nR/(sum(Vi/Ti))
    \f]    
    where the i might stand for the gap, plena, central hole or cracks.
    The calculation of the various V/T terms is done separately for the various
    section of the free volume. 

    The calculation of the gap and plena volumes is done applying the 
    gauss-green theorem applied to the bounding surfaces:
    \f[
        volume = volIntegral(dV) = 1/3*volIntegral(div(r) dV) = 1/3*surfIntegral(S*n dS)
    \f]    

    Usually, the cladding inner surface is larger than the fuel outer surface.
    Also, the mesh is not necessarily conformal on the two sides of the gap.
    Therefore it is necessary to introduce a scaling factor in order to 
    separate the plena from the gap volume. 

    The scaling factor on a owner face is equal to 1 if the face is entirely 
    covered by the opposing surface; It is 0 if it does not see any portion of 
    the opposing surface; It is a fraction of 1 if it is only partly covered by 
    the nbr surface. By definition, the correction factor for the fuel should
    always be equal to 1. For the cladding they should be equal to 1 for the 
    central section facing the fuel; 0 for the faces facing the plena; a 
    fractional value for those faces that are partially facing the fuel and
    partially facing the plenum.

    The correction factors could be derived from the weightSum function of 
    the AMI, which by definition provides the ratio between the amount of src 
    surface projected on the tgt face face divided by the src/tgt face.
    However, because the rod is a cylindrical body, the weightSum does not 
    differentiate between two faces separated by a gap (i.e. in a cylinder they 
    do not have the same area) and two faces that are partially overlapping on 
    the axial direction.

    For this reason, we find the scaling factor in an alternative manner, using
    cutting planes and costructing the overlapping faces manually with the 
    points found by intersecting the cutting planes and the edges of cladding 
    inner patch faces.

Usage
    In solverDict file:
    \verbatim
    gapGas FRAPCON;

    gapGasOptions
    {
        gapPatches ( fuelOuter cladInner );
        holePatches ( fuelInner );
        topFuelPatches ( fuelTop );
        bottomFuelPatches ( fuelBottom );
        
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

        includeCentralHole on;
        includeDishes on;
    }
    \endverbatim

    In gapGas file inside the time step folder (e.g. 0/):    
    \verbatim
    
    massFractions
    {
        Ar 0;
        He 1;
        Kr 0;
        Ne 0;
        Rn 0;
        Xe 0;
    };

    gasPressureType                 fromModel;// fixed;// fromList;
    gasPressure                     2.25e6;
     
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
    gapFRAPCON.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\contribution
    E. Brunetto, C. Fiorina - EPFL

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef gapFRAPCON_H
#define gapFRAPCON_H

#include "gapGasModel.H"
#include "Function1.H"
#include "globalOptions.H"
#include <set>
#include "Table.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class gapFRAPCON Declaration
\*---------------------------------------------------------------------------*/

class gapFRAPCON
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

        //- Correct scaling factors
        virtual void correctScalingFactors();

        //- Calculate average axial location from a list of patches
        virtual scalar patchListAvgAxialLocation(const std::set<label>&);

        //- Estimate the current gap quantities (volume, temperature, VoverT)
        virtual void correctGap();  

        //- Estimate the current hole quantities (volume, temperature, VoverT)
        virtual void correctHole();  

        //- Calculate the plena V, T and V/T
        virtual void correctPlena();
        
        //- Calculate the T, V, and the sum of V/T for gas between dishes
        virtual void correctDish();
        
        //- Calculate the crack T, V, and the sum of V/T for gas in cracks in fuel
        virtual void correctCrack();
        
        //- Disallow default bitwise copy construct
        gapFRAPCON(const gapFRAPCON&);

        //- Disallow default bitwise assignment
        void operator=(const gapFRAPCON&);
        
protected:  

    // Protected data 

        //- Dish volume [m3]
        scalar dishV_;    

        //- Hole volume [m3]
        scalar holeV_;    

        //- Top plenum volume [m3]
        scalar topPlenumV_;   

        //- Bottom plenum volume [m3]
        scalar bottomPlenumV_;   

        //- Cracks volume [m3]
        scalar crackV_;   

        //- Gap volume offset [m3]   
        scalar gapOffsetV_;  

        //- Additional gas reserve volume [m3]
        scalar reserveV_;

        //- Dish average temperature [K]
        scalar dishT_;

        //- Hole average temperature [K]
        scalar holeT_;

        //- Top plenum average temperature [K]
        scalar topPlenumT_;

        //- Bottom plenum average temperature [K]
        scalar bottomPlenumT_;

        //- Crack gas average temperature [K]
        scalar crackT_;

        //- Additional gas reserve temperature [K]   
        scalar reserveT_; 

#ifdef OPENFOAMFOUNDATION
        //- Additional gas reserve temperature (time list) [K]
        autoPtr<Function1s::Table<scalar>> reserveTlist_;  
#elif OPENFOAMESI        
        //- Additional gas reserve temperature (time list) [K]
        autoPtr<Function1<scalar>> reserveTlist_;   
#endif 

        //- Gap sum(Vi/Ti) for each face of the gap [m3/K]
        scalar gapVoverT_;

        //- Sum(Vi/Ti) for the dish volume [m3/K]
        scalar dishVoverT_;

        //- Hole sum(Vi/Ti) for each face of the hole [m3/K]
        scalar holeVoverT_;

        //- Crack sum(Vi/Ti) for each fuel cell [m3/K]
        scalar crackVoverT_;

        //- Highest axial coordinate of fuel region
        scalar zMax_;
        
        //- Lowest axial coordinate of fuel region
        scalar zMin_;
        
        //- List of scaling factors for the gap (one list per patch)
        List<scalarField> scalingFactors_;
        
        //- List of scaling factors for top plenum (one list per patch)
        List<scalarField> scalingFactorsTopPlenum_;
        
        //- List of scaling factors for bottom plenum (one list per patch)
        List<scalarField> scalingFactorsBottomPlenum_;

        //- Include fuel central hole volume or not
        Switch includeCentralHole_;

        //- Hash table of patch groups bounding the gap (gap, fuelTop, fuel bottom)
        //  TODO: separate fuelGap and claddingGap
        autoPtr<Foam::HashTable<std::set<label>, Foam::word>> patchIDsPtr_;

        //- Name of relocation strain field
        word relocationName_;
        
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

        //- Axial pin direction
        vector pinDirection_;

public:           

    //- Runtime type information
    TypeName("FRAPCON");
    

    // Constructors

        //- Construct from mesh, materials and dict
        gapFRAPCON
        (   
            const fvMesh& mesh,
            const materials& mat,
            const volScalarField& T,
            const volVectorField& DorDD,
            const fissionGasRelease& fgr,
            const dictionary& gapGasOptDict
        );
        

    //- Destructor
    ~gapFRAPCON();


    // Member Functions
            
        //- Gap gas thermal conductivity (from FRAPCON manual)
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
