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
    Foam::phaseTransitionZyDynamic

Description
    Class for the beta phase transition of Zy-based claddings, dynamic approach
    derived from "A.R. Massih, Transformation kinetics of zirconium alloys under
    non-isothermal conditions, Journal of Nuclear Materials, Volume 384, Issue 3,
    2009, Pages 330-335, ISSN 0022-3115, https://doi.org/10.1016/j.jnucmat.2008.11.033.
    (https://www.sciencedirect.com/science/article/pii/S0022311508007654).

Usage
    In solverDict file:
    \verbatim
    materials
    {
        ...

        cladding
        {
            material zircaloy;

            ...

            phaseTransitionModel     ZyDynamic;
        }

        ...
    }
    \endverbatim

SourceFiles
    phaseTransitionZyDynamic.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef phaseTransitionZyDynamic_H
#define phaseTransitionZyDynamic_H

#include "phaseTransitionModel.H"
#include "physicoChemicalConstants.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class phaseTransitionZyDynamic Declaration
\*---------------------------------------------------------------------------*/

class phaseTransitionZyDynamic
:
    public phaseTransitionModel
{
    // Private data

        //- Hydrogen concentration 
        // TODO (EB) : this should be provided by a dedicated model.   
        const scalar w_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        phaseTransitionZyDynamic(const phaseTransitionZyDynamic&);

        //- Disallow default bitwise assignment
        void operator=(const phaseTransitionZyDynamic&);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("ZyDynamic");


    // Declare run-time constructor selection table


    // Constructors

        //- Construct from dictionary
        phaseTransitionZyDynamic
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Selectors


    //- Destructor
    virtual ~phaseTransitionZyDynamic();


    // Member Functions
    
        //- Update beta fraction
        virtual void correct(const scalarField& T, const labelList& addr);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //