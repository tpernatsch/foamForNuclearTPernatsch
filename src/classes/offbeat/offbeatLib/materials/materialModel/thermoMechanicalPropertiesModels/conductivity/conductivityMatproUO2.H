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
    Foam::conductivityMatproUO2

Description
    Model for UO2 conductivity from MATPROv11.

SourceFiles
    conductivityMatproUO2.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef conductivityMatproUO2_H
#define conductivityMatproUO2_H

#include "conductivityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class conductivityMatproUO2 Declaration
\*---------------------------------------------------------------------------*/

class conductivityMatproUO2
:
    public conductivityModel
{
    // Private data
        
    //- Name used for the burnup field
    const word burnupName_;
    
    //- Reference to Burnup field
    const volScalarField* Bu_;

    //- Fraction of Theoretical Density 
    scalar densityFrac_;

    //- Initial content of Gadolinium
    scalar GdContent_;

    //- Parameters for the Matpro model
    scalar par1;
    scalar par2;
    scalar par3;
    scalar par4;
    scalar par5;
    scalar par6;
    scalar par7;
    scalar par8;
    scalar par9;
    scalar par10;
    scalar par11;
    scalar par12;
    scalar par13;
    
    //- Perturbation factor
    scalar perturb;

    // Private Member Functions

    //- Disallow default bitwise copy construct
    conductivityMatproUO2(const conductivityMatproUO2&);

    //- Disallow default bitwise assignment
    void operator=(const conductivityMatproUO2&);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("UO2MATPRO");

    // Constructors

    //- Construct from dictionary
    conductivityMatproUO2
    (
        const fvMesh& mesh,
        const dictionary& dict, 
        const word defaultModel   
    );
    
    //- Destructor
    virtual ~conductivityMatproUO2();


    // Member Functions
    
    //- Update conductivity  
    virtual void correct(scalarField& sf, const scalarField& T, const labelList& addr);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
