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
    Foam::swellingFBRMOX

Description
    Class modelling swelling strain for FBR MOX fuel from Dienst at al.
    Sources : 
    - Dienst, W., Muelle-Lyda, I., Zimmerman, H. Swelling, densification and 
    creep of oxide and carbide fuels under irradiation, in: International 
    conference on fast breeder reactor performance, 5-8 March 1979, Monterey, 
    California, USA, 1979, pp. 166–175. 
    - "Validazione e analisi dei gas di fissione in combustibili MOX ad 
    elevato burnup", D.Rozzia, N.Forgione, A. Ardizzone, A. Del Nevo.

Usage
    In solverDict file:
    \verbatim
    materials
    {
        fuel
        {
            material UPuO2;

            ...

            swellingModel     UPuO2Fbr;
        }

        ...
    }
    \endverbatim

SourceFiles
    swellingFBRMOX.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro - EPFL\n

\date
    July 2024
\*---------------------------------------------------------------------------*/

#ifndef swellingFBRMOX_H
#define swellingFBRMOX_H

#include "swellingModel.H"
#include "sliceMapper.H"
#include "physicoChemicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class swellingFBRMOX Declaration
\*---------------------------------------------------------------------------*/

class swellingFBRMOX
:
    public swellingModel
{
    // Private data

        //- Name used for the burnup field
        const word burnupName_;

        //- Reference to Burnup field
        const volScalarField* Bu_;

        //- Name used for the gapWidth field
        const word gapWidthName_;

        //- Reference to gapWidth field
        const volScalarField* gapWidth_;

        //- Fuel outersurface patchID to estimate average gapwidth
        label outerPatchID_;

        //- Parameters for the Frapcon model
        scalar par1;
        scalar par2;
        scalar par3;

        //- Reference to slice mapper
        const sliceMapper* mapper_;

        // -
        const volSymmTensorField* intraGasSwell_;
        const volSymmTensorField* interGasSwell_;

    // Private Member Functions

        //- Compute gapWidth 1D to assess if gap is closed for each slice
        virtual scalarField computeGapWidth1D();

        //- Disallow default bitwise copy construct
        swellingFBRMOX(const swellingFBRMOX&);

        //- Disallow default bitwise assignment
        void operator=(const swellingFBRMOX&);

    // Scaling parameters for UQ and SA
        
        //- Group name for parameter definition
        declareNameWithVirtual(swellingFBRMOX, group_, group);
    
        //- Swelling strain scaling factor [-]
        declareParameterWithVirtual(swellingFBRMOX, 
            F_epsilonSwelling_, F_epsilonSwelling);
        
        //- Swelling strain offset [-]
        declareParameterWithVirtual(swellingFBRMOX, 
            delta_epsilonSwelling_, delta_epsilonSwelling);

protected:

    // Protected data

public:

    //- Runtime type information
    TypeName("UPuO2Fbr");

    // Declare run-time constructor selection table

    // Constructors

        //- Construct from dictionary
        swellingFBRMOX
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dictionary& baseDict
        );


    // Selectors

    //- Destructor
    virtual ~swellingFBRMOX();


    // Member Functions

    //- Update heatCapacity
    virtual void correct(const scalarField& T, const labelList& addr) ;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
