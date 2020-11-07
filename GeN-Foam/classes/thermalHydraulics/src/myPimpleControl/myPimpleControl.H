/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Released 2011 OpenCFD Ltd.
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Modified code Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::myPimpleControl

Description
    PIMPLE control class to supply convergence information/checks for
    the PIMPLE loop. In addition to standard PIMPLE functions, this class
    adds support for a variable number of PISO iterations (via the 
    corrPISOUntilConvergence_ flag) within each PIMPLE iteration, and
    adds an optional residual control for the number of PIMPLE iterations
    performed in each time step based on the evaulation of the initialResidual
    of the first PISO iteration in every PIMPLE iteration (rather than the
    latest available initialResidual, from any PISO iteration, as done by the
    base class). This proves particularly more useful for multiphase
    simulation.

Author
    Stefan Radman (stefan.radman@epfl.ch; stefanradman92@gmail.com)

Dictionary example (specifically to show new entries use)

PIMPLE
{
    nOuterCorrectors            12;
    nCorrectors                 8;
    nNonOrthogonalCorrectors    0;
    corrPISOUntilConvergence    true;   //- (new entry, default: false)

    residualControl
    {
        p_rgh
        {
            tolerance                   5e-4;
            relTol                      0;
            useFirstPISOInitialResidual true;   //- (new entry, default: false)
        }
    }
}

\*---------------------------------------------------------------------------*/

#ifndef myPimpleControl_H
#define myPimpleControl_H

#include "pimpleControl.H"

//- Declare that myPimpleControl will be used
#define MYPIMPLE_CONTROL

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class myPimpleControl Declaration
\*---------------------------------------------------------------------------*/

class myPimpleControl
:
    public pimpleControl
{
    // Private member functions

        //- No copy construct
        myPimpleControl(const myPimpleControl&) = delete;

        //- No copy assignment
        void operator=(const myPimpleControl&) = delete;


protected:

    // Protected data

        // Solution controls

            //- Minimum amount of PIMPLE iterations to be performed
            int minNCorrPIMPLE_;

            //- If true, the number of PISO iterations performed in every
            //  PIMPLE iteration will be at most nCorrPISO_, but the PISO
            //  iterations will stop if the linear solver has reached 
            //  convergence. This operates completely independently form the
            //  novel optional residualControl for pressure based on 
            //  first PISO initial residuals
            bool corrPISOUntilConvergence_;

            //-
            int nCorrPISOInPrevPIMPLE_;

            //-
            int nCorrPISOInPrevPrevPIMPLE_;

            //-
            bool stopLoop_;

    // Protected Member Functions

        //- Based on maxTypeResidual in solutionControl. The idea is to get 
        //  the initResidual of the first PISO corrector loop of the previous
        //  PIMPLE iteration. This is because it is a good indicator of the
        //  convergence of pressure-velocity coupling in multiphase 
        //  scenarios. This and the next function work in tandem. In
        //  particular, the residual pair consists of the initial residual of
        //  the very first PISO corrector of the first PIMPLE iteration, while
        //  the second contains the initial residual of the first PISO
        //  corrector of the previous PIMPLE iteration
        template<class Type>
        bool firstPISOPrevPIMPLETypeResidual
        (
            const entry& solverPerfDictEntry,
            Pair<scalar>& residual
        ) const;

        //-
        Pair<scalar> firstPISOPrevPIMPLEResidual
        (
            const entry& solverPerfDictEntry
        ) const;

        //- Read controls from fvSolution dictionary
        virtual bool read();

        //- Return true if all convergence checks are satisfied. This modified
        //  version can used the first PISO initialResidual for the evaulation
        //  of convergence criteria if the useFirstPISOInitialResidualFlag
        //  is specified in the field residualControl dict. 
        //  It is meant specifically for the pressure field and does not make
        //  sense to use it with any other field. Nonetheless, the flag 
        //  defaults to false if not provided and this function behaves exactly
        //  like its original version
        virtual bool criteriaSatisfied();

public:

    // Static Data Members

        //- Run-time type information
        TypeName("myPimpleControl");


    // Constructors

        //- Construct from mesh and the name of control sub-dictionary
        myPimpleControl
        (
            fvMesh& mesh,
            const word& dictName = "PIMPLE",
            const bool verbose = true
        );


    //- Destructor
    virtual ~myPimpleControl() = default;


    // Member Functions

        // Access

        // Solution control

            //- PIMPLE loop
            virtual bool loop();

            //- Pressure corrector loop control
            inline bool correct();

            //- Pressure corrector loop control with interrupt control
            inline bool correct(bool interrupt);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "myPimpleControlI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
