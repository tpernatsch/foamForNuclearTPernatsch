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
    Wrapping class for an fgr model coupling OFFBEAT to SCIANTIX.

    \htmlinclude fissionGasRelease_fgrSCIANTIX.html

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


    // Private Member Functions

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
          const dictionary& fgrDict
        );


    //- Destructor
    ~fgrSCIANTIX();


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