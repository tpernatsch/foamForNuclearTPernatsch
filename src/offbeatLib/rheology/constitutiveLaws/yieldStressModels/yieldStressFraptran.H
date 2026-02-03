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
    Foam::yieldStressFraptran

Description
    yield stress model from WGFS RIA fuel codes benchmark Phase II specifications. 
    https://www.oecd-nea.org/jcms/pl_19704/reactivity-initiated-accident-ria-fuel-codes-benchmark-phase-ii-volume-2-task-no-1-specifications?details=true

    The model itself comes from the Fraptran fuel code and the Matpro database.
    The yield stress is dependant on temperature, fast neutronce fluence, Cold Work amount and strain rate.

    \htmlinclude yieldStressModel_yieldStressFraptran.html

SourceFiles
    yieldStressFraptran.C

\mainauthor
    Matthieu Reymond

\contribution

\date 
    Feb. 2023


\*---------------------------------------------------------------------------*/

#ifndef yieldStressFraptran_H
#define yieldStressFraptran_H

#include "yieldStressModel.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "zeroGradientFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class "yieldStressConstant" Declaration
\*---------------------------------------------------------------------------*/

class yieldStressFraptran
:
    public yieldStressModel 
{    

        // Private data    
        // Private Member Functions
        //- Disallow default bitwise copy construct
        yieldStressFraptran(const yieldStressFraptran&);

        //- Disallow default bitwise assignment
        void operator=(const yieldStressFraptran&);

        scalar phiValue_;

        scalar CWValue_;
        
protected:
        
        // Protected data
        volScalarField& Kcoeff_;

        volScalarField& m_;

        volScalarField& n_;

//        volSymmTensorField EpsDot_;
    
public:
 
    //- Runtime type information
    TypeName("Fraptran");

    // Constructors

        //- Construct from dictionary
        yieldStressFraptran
        (
            const fvMesh& mesh,
            const dictionary& lawDdict
        );

    //- Destructor
   virtual ~yieldStressFraptran();

    // Member Functions
    
        //- Update yield stress
        virtual void correctYieldStress
        (
            const labelList& addr
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
